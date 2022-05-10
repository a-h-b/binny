"""
Created on Wed Feb 22 10:50:35 2021

@author: oskar.hickl
"""
import itertools
import logging
from pathlib import Path
from timeit import default_timer as timer
import re

import hdbscan
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import parallel_backend, Parallel, delayed
from mpl_toolkits.mplot3d import Axes3D
from openTSNE import TSNEEmbedding, affinity
from openTSNE import initialization
from skbio.stats.composition import clr, multiplicative_replacement
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier

def unify_multi_model_genes(gene, markers='essential'):
    # Dict to unify genes with multiple models
    if markers == 'essential':
        hmm_dict = {'TIGR00388': 'glyS', 'TIGR00389': 'glyS', 'TIGR00471': 'pheT', 'TIGR00472': 'pheT',
                    'TIGR00408': 'proS', 'TIGR00409': 'proS', 'TIGR02386': 'rpoC', 'TIGR02387': 'rpoC'}
        # Replace gene with unified name if it has multiple models
        if gene in hmm_dict:
            gene = hmm_dict.get(gene)
    return gene


def gff2ess_gene_df(annotation_file, target_attribute='essential', get_dict=False):
    contig_dict = {}
    # Get all info for each line
    with open(annotation_file, 'r') as af:
        for line in af:
            line = line.strip(' \t\n').split('\t')
            contig, contig_start, contig_end, attributes = line[0], line[3], line[4], line[8].split(';')
            # Check if feature has no start and/or stop position or length = 0
            if not contig_start or not contig_end or int(contig_start) - int(contig_end) == 0:
                continue
            # Only add line if feature was annotated as essential gene
            for attribute in attributes:
                attribute_name, attribute_id = attribute.split('=')[0], attribute.split('=')[1]
                if attribute_name == target_attribute and attribute_id:
                    if not contig_dict.get(contig):
                        if target_attribute == 'essential':
                            contig_dict[contig] = [unify_multi_model_genes(attribute_id)]
                        elif target_attribute == 'checkm_marker':
                            contig_dict[contig] = [marker.split('.')[0] for marker in attribute_id.split(',')]
                        else:
                            contig_dict[contig] = [attribute_id]
                    else:
                        if target_attribute == 'essential':
                            contig_dict[contig].append(unify_multi_model_genes(attribute_id))
                        elif target_attribute == 'checkm_marker':
                            contig_dict[contig].extend([marker.split('.')[0] for marker in attribute_id.split(',')])
                        else:
                            contig_dict[contig].append(attribute_id)
    # Create list of lists to join dict items as string and build data frame
    annot_cont_line_list = [[contig, ','.join(contig_dict.get(contig))] for contig in contig_dict]
    annot_cont_df = pd.DataFrame(annot_cont_line_list, columns=['contig', 'essential'], dtype='string[pyarrow]')
    if get_dict:
        return annot_cont_df, contig_dict
    else:
        return annot_cont_df


def load_and_merge_cont_data(annot_file, depth_file, coords, dims, coords_from_file=True):
    logging.info('Reading annotation.')
    annot_df = gff2ess_gene_df(annot_file, target_attribute='checkm_marker')
    logging.info('Reading embedding coordinates.')
    dim_range = [i + 1 for i in range(dims)]
    if coords_from_file:
        coord_cols = ['contig'] + ['dim' + str(i) for i in dim_range]
        coord_dtypes = ['string[pyarrow]'] + ['float32' for i in dim_range]
        coord_dtype_dict = {col: dtype for col, dtype in zip(coord_cols, coord_dtypes)}
        coord_df = pd.read_csv(coords, sep='\t', names=coord_cols, dtype=coord_dtype_dict)
    else:
        coord_df = coords

    with open(depth_file, 'r') as f:
        for line in f:
            first_line = line.strip(' \t\n').split('\t')
            n_depths = len(first_line[1:])
            break
    logging.info(f'{n_depths} samples for depth data found.')

    depth_cols = ['contig'] + ['depth' + str(depth_ind + 1) for depth_ind in range(n_depths)]
    depth_dtypes = ['string[pyarrow]'] + ['float32' for i in range(n_depths)]
    depth_dtype_dict = {col: dtype for col, dtype in zip(depth_cols, depth_dtypes)}
    depth_df = pd.read_csv(depth_file, sep='\t', low_memory=False, names=depth_cols, dtype=depth_dtype_dict)

    annot_df['contig'] = annot_df['contig'].astype('string')
    coord_df['contig'] = coord_df['contig'].astype('string')
    depth_df['contig'] = depth_df['contig'].astype('string')

    # Merge annotation data and coords_file first, keeping only contigs with coords, then merge with depth data
    cont_data_df = annot_df.merge(coord_df, how='right', on='contig') \
        .merge(depth_df, how='inner', on='contig').sort_values(by='contig')
    downcast_pd_df(cont_data_df)

    return cont_data_df


def knn_sne_coords(contig_info_df, knn_sne_pk, dims):
    if knn_sne_pk < 2 * len(dims):
        knn_sne_pk = 2 * len(dims)
    if knn_sne_pk > contig_info_df['contig'].size:
        knn_sne_pk = contig_info_df['contig'].size - 1
    if knn_sne_pk <= 0:
        knn_sne_pk = 1

    knn_c = KNeighborsClassifier(n_neighbors=knn_sne_pk, weights='distance')
    coords_df = contig_info_df.loc[:, dims]

    knn_c.fit(coords_df, contig_info_df.loc[:, 'contig'])
    # sort distance of neighbouring points
    try:
        bin_kneighbors = knn_c.kneighbors(n_neighbors=knn_sne_pk)
    except ValueError:
        logging.warning('knn_c', knn_c)
        logging.warning('knn_sne_pk', knn_sne_pk)
        logging.warning('contig_info_df', contig_info_df)
        logging.warning('contig_info_df.size', contig_info_df.size)
        logging.warning('contig_info_df[\'contig\']', contig_info_df['contig'])
        logging.warning('contig_info_df[\'contig\'].size', contig_info_df['contig'].size)
        logging.warning('coords_df', coords_df)
        logging.warning('coords_df.size', coords_df.size)
        if knn_sne_pk > 1:
            bin_kneighbors = knn_c.kneighbors(n_neighbors=knn_sne_pk - 1)
        else:
            bin_kneighbors = knn_c.kneighbors(n_neighbors=1)
    skn = pd.Series(sorted([row[-1] for row in bin_kneighbors[0]]), dtype='float32')
    # calculate running standard deviation between 10 neighbouring points
    sdkn = skn.rolling(10, center=True, min_periods=0).std()
    # find the first jump in distances at the higher end
    try:
        est = sorted([e for i, e in zip(sdkn, skn) if i > sdkn.quantile(.975) and e > skn.median()])[0]
    except IndexError:
        est = None
    return est


def contig_df2cluster_dict(contig_info_df, dbscan_labels, use_noise=False):
    # This function is garbage. From pandas df to dict to df. Need to improve asap
    if use_noise:
        use_noise = 'noise_stays'
    else:
        use_noise = '-1'
    cluster_dict = {}

    dims = [i for i in contig_info_df.columns if 'dim' in i]
    depths = [i for i in contig_info_df.columns if 'depth' in i]
    contig_info_df = contig_info_df.loc[:, ['contig', 'essential'] + depths + dims]
    contig_info_data = ['contigs', 'essential'] + depths + dims
    contig_info_df['cluster'] = dbscan_labels
    contig_info_df['cluster'] = contig_info_df['cluster'].astype('string[pyarrow]')

    contig_info_df['essential'] = contig_info_df['essential'].astype('string').fillna('non_essential').astype('string[pyarrow]')

    try:
        tmp_contig_dict = contig_info_df.set_index('contig').to_dict('index')
    except ValueError:
        print("contig_df2cluster_dict, contig_info_df['contig']:", len(contig_info_df['contig'].to_list()))
        print("contig_df2cluster_dict, contig_info_df['contig']:", len(set(contig_info_df['contig'].to_list())))
        raise Exception

    for contig, contig_data in tmp_contig_dict.items():
        # contig_cluster = shorten_cluster_names(contig_data.get('cluster'))
        contig_cluster = shorten_cluster_names(contig_data['cluster'])
        contig_data_list = [contig] + [contig_data.get(i) for i in contig_info_data[1:]]
        if not cluster_dict.get(contig_cluster) and use_noise not in contig_cluster:
            cluster_dict[contig_cluster] = {info: [data] for info, data in zip(contig_info_data, contig_data_list)}
        elif use_noise not in contig_cluster:
            for info, data in zip(contig_info_data, contig_data_list):
                cluster_dict.get(contig_cluster, {}).get(info).append(data)

    return cluster_dict


def hdbscan_cluster(contig_data_df, pk=None, include_depth=False, n_jobs=1, hdbscan_epsilon=0.25, hdbscan_min_samples=2,
                    dist_metric='manhattan'):
    dims = [i for i in contig_data_df.columns if 'dim' in i and '_' not in i]
    depths = [i for i in contig_data_df.columns if 'depth' in i]

    if not include_depth:
        dim_df = contig_data_df.loc[:, dims]
    else:
        dim_df = contig_data_df.loc[:, dims + depths]
    if not pk:
        pk = get_perp(contig_data_df['contig'].size)
        if pk < len(dims) * 2:
            pk = len(dims) * 2

    with parallel_backend('threading'):
        np.random.seed(0)
        logging.debug(f'HDBSCAN params: min_cluster_size={pk}, min_samples={hdbscan_min_samples}, '
                     f'cluster_selection_epsilon={hdbscan_epsilon}, metric={dist_metric}.')

        hdbsc = hdbscan.HDBSCAN(core_dist_n_jobs=n_jobs, min_cluster_size=pk, min_samples=hdbscan_min_samples,
                                cluster_selection_epsilon=hdbscan_epsilon, metric=dist_metric).fit(dim_df)

    cluster_labels = hdbsc.labels_

    while len(set(cluster_labels)) == 1 and str(list(set(cluster_labels))[0]) == '-1' and pk >= len(dims) * 2:
        pk = int(pk * 0.75)
        logging.warning('HDBSCAN found only noise trying with lower min_cluster_size={0}.'.format(pk))
        with parallel_backend('threading'):
            np.random.seed(0)
            hdbsc = hdbscan.HDBSCAN(core_dist_n_jobs=n_jobs, min_cluster_size=pk, min_samples=hdbscan_min_samples,
                                    cluster_selection_epsilon=hdbscan_epsilon, metric=dist_metric).fit(dim_df)
        cluster_labels = hdbsc.labels_

    if len(set(cluster_labels)) == 1 and str(list(set(cluster_labels))[0]) == '-1':
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels, use_noise=True)
    else:
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)

    return cluster_dict, cluster_labels


def run_initial_scan(contig_data_df, initial_cluster_mode, dbscan_threads, pk=None, include_depth=False,
                     hdbscan_epsilon=0.25, hdbscan_min_samples=2, dist_metric='manhattan'):

    if not pk:
        pk = get_perp(contig_data_df['contig'].size)
    first_clust_dict, labels = {}, []
    if initial_cluster_mode == 'HDBSCAN' or not initial_cluster_mode:
        first_clust_dict, labels = hdbscan_cluster(contig_data_df, pk=pk, n_jobs=dbscan_threads,
                                                   include_depth=include_depth, hdbscan_epsilon=hdbscan_epsilon,
                                                   hdbscan_min_samples=hdbscan_min_samples, dist_metric=dist_metric)
    return first_clust_dict, labels


def cluster_df_from_dict(cluster_dict):
    cluster_df = pd.DataFrame()
    if len(set(list(cluster_dict.keys()))) == 1 and str(list(set(list(cluster_dict.keys())))[0]) == '-1':
        cluster_df['cluster'] = [cluster for cluster in cluster_dict]
    else:
        cluster_df['cluster'] = pd.Series([cluster for cluster in cluster_dict if not cluster == '-1'],
                                          dtype='string[pyarrow]')
    for metric in ['contigs', 'essential', 'completeness', 'purity']:
        metric_list = []
        metric_uniq_list = None
        if metric == 'essential':
            metric_uniq_list = []
        for cluster in cluster_dict:
            if metric in ['contigs', 'essential']:
                data_list = [i for e in cluster_dict.get(cluster, {}).get(metric) for i in e.split(',')
                             if not e == 'non_essential']
                metric_list.append(len(data_list))
                if metric == 'essential':
                    metric_uniq_list.append(len(set(data_list)))
            else:
                metric_list.append(cluster_dict.get(cluster, {}).get(metric, 0))
        if metric in ['contigs', 'essential']:
            cluster_df[metric] = pd.Series(metric_list, dtype='string[pyarrow]')
        else:
            cluster_df[metric] = pd.Series(metric_list, dtype='float32')
        if metric_uniq_list:
            try:
                cluster_df['unique_' + metric] = pd.Series(metric_uniq_list, dtype='string[pyarrow]')
            except ValueError:
                logging.warning(metric_uniq_list, metric, cluster_df, cluster_dict)
                raise Exception
    return cluster_df.fillna(0)


def cluster_dict2bin_dict(cluster_dict, cluster_df):
    cluster_df = cluster_df.set_index('cluster')
    bin_dict = {}
    for cluster in cluster_dict:
        cluster_purity = str(int(cluster_df.loc[cluster, 'purity'] * 100))
        cluster_completeness = str(int(cluster_df.loc[cluster, 'completeness'] * 100))
        bin_name = '_'.join(['binny'] + [cluster] + [cluster_purity] + [cluster_completeness])
        bin_dict[bin_name] = cluster_dict.get(cluster, {}).get('contigs')
    return bin_dict


def sort_cluster_dict_data_by_depth(cluster_dict):
    sorted_dict = {}
    dims = [i for i in list(cluster_dict[list(cluster_dict.keys())[0]].keys()) if 'dim' in i]
    depths = [i for i in list(cluster_dict[list(cluster_dict.keys())[0]].keys()) if 'depth' in i]
    for cluster in cluster_dict:
        sorted_dict[cluster] = {}
        deps = np.array(cluster_dict[cluster]['depth1']).argsort()
        for metric in depths + ['contigs', 'essential'] + dims:
            metric_np = np.array(cluster_dict[cluster][metric])
            sorted_dict[cluster][metric] = metric_np[deps]
    return sorted_dict


def gather_cluster_data(cluster, cluster_dict, marker_sets_graph, tigrfam2pfam_data_dict):
    cluster_essential_genes = [gene for genes in cluster_dict.get(cluster, {}).get('essential')
                               for gene in genes.split(',') if not gene == 'non_essential']
    if cluster_essential_genes:
        marker_set = choose_checkm_marker_set(cluster_essential_genes, marker_sets_graph, tigrfam2pfam_data_dict, max_depth_lvl=4)
        taxon, cluster_completeness, cluster_purity = marker_set[0], round(marker_set[1], 3), round(marker_set[2], 3)
    else:
        cluster_purity = 0
        cluster_completeness = 0
        taxon = 'none'
    cluster_info = sorted(list(cluster_dict.get(cluster).keys()))
    cluster_data = [cluster_dict.get(cluster, {}).get(i) for i in cluster_info] + [cluster_purity, cluster_completeness,
                                                                                   taxon]
    return cluster_data


def hdbscan_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan, depth=False, max_tries=100,
                         hdbscan_epsilon=0.25, hdbscan_min_samples=2, dist_metric='manhattan'):
    start = timer()

    dbscan_tries = 0
    new_clusters_labels = [1]

    logging.debug('Working on {0}.'.format(shorten_cluster_names(cluster)))

    while len(set(new_clusters_labels)) == 1 and dbscan_tries <= max_tries:
        dbscan_tries += 1

        dims = [i for i in cluster_contig_df.columns if 'dim' in i and '_' not in i]
        if depth:
            depths = [i for i in cluster_contig_df.columns if 'depth' in i]
            dims += depths

        df_coords = cluster_contig_df.loc[:, dims].to_numpy(dtype=np.float32)

        if cluster_contig_df['contig'].size == 1:
            new_clusters = contig_df2cluster_dict(cluster_contig_df, [cluster + '.' + str(0)])
            logging.debug(
                'Cluster {0} is only a single data point. Returning it.'.format(shorten_cluster_names(cluster)))
            return [new_clusters, cluster]

        with parallel_backend('threading'):
            np.random.seed(0)
            hdbsc = hdbscan.HDBSCAN(core_dist_n_jobs=threads_for_dbscan, min_cluster_size=pk,
                                    min_samples=hdbscan_min_samples, cluster_selection_epsilon=hdbscan_epsilon,
                                    metric=dist_metric).fit(df_coords)
        new_clusters_labels = hdbsc.labels_
        if len(set(new_clusters_labels)) > 1:
            new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
                                 enumerate(set(new_clusters_labels))}
            new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
            new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
            end = timer()
            logging.debug('Found {0} sub-clusters with HDBSCAN in cluster {1} after {2}s'.format(
                len(set(new_clusters_labels)), shorten_cluster_names(cluster), int(end - start)))
            return [new_clusters, cluster]

        else:
            end = timer()
            logging.debug('Failed to find sub-clusters for {0} with HDBSCAN after {1}s.'.format(
                shorten_cluster_names(cluster), int(end - start)))
            return [{}, cluster]
    end = timer()
    logging.debug(
        'Failed to find sub-clusters for {0} with HDBSCAN after {1}s.'.format(
            shorten_cluster_names(cluster), int(end - start)))
    return [{}, cluster]


def get_sub_clusters(cluster_dicts, threads_for_dbscan, marker_sets_graph, tigrfam2pfam_data_dict, purity_threshold=0.95,
                     completeness_threshold=0.9, pk=None, cluster_mode=None, include_depth=True, hdbscan_epsilon=0.25,
                     hdbscan_min_samples=2, dist_metric='manhattan'):
    outputs = []

    for cluster_dict in cluster_dicts:
        cluster = list(cluster_dict.keys())[0]
        # Create dict with just cluster and sort again, to ensure order by depth is
        cluster_dict = sort_cluster_dict_data_by_depth({cluster: cluster_dict[cluster]})
        # All data needed stored in list with following order:
        # Indices for clust_dat 0-1 are always:
        # contigs and depth, then n dimensions, -4 to -1 are always essential, purity, completeness, selected taxon
        clust_dat = gather_cluster_data(cluster, cluster_dict, marker_sets_graph, tigrfam2pfam_data_dict)
        clust_pur = float(clust_dat[-3])
        clust_comp = float(clust_dat[-2])
        clust_taxon = clust_dat[-1]

        cluster_dict[cluster]['purity'] = clust_pur
        cluster_dict[cluster]['completeness'] = clust_comp
        cluster_dict[cluster]['taxon'] = clust_taxon

        cluster_pur_thresh = purity_threshold

        if 0.700 < clust_comp <= 0.800:
            if cluster_pur_thresh < 0.925:
                cluster_pur_thresh = 0.925
        elif clust_comp <= 0.700:
            if cluster_pur_thresh < 0.95:
                cluster_pur_thresh = 0.95

        if clust_taxon in ['Bacteria', 'Archaea']:
            if 0.850 < clust_comp <= 0.900:
                if cluster_pur_thresh < 0.925:
                    cluster_pur_thresh = 0.925
            elif 0.800 < clust_comp <= 0.850:
                if cluster_pur_thresh < 0.95:
                    cluster_pur_thresh = 0.95
            elif 0.700 < clust_comp <= 0.800:
                if cluster_pur_thresh < 0.975:
                    cluster_pur_thresh = 0.975
            elif clust_comp <= 0.700:
                if cluster_pur_thresh < 0.99:
                    cluster_pur_thresh = 0.99

        if clust_pur < cluster_pur_thresh and isinstance(clust_pur, float) and clust_comp >= completeness_threshold:
            logging.debug('Cluster {0} below purity of {1} with {2} and matches completeness of {3} with {4}. '
                          'Attempting to split.'.format(shorten_cluster_names(cluster), purity_threshold, clust_pur,
                                                        completeness_threshold, clust_comp))

            # initialise some stuff
            cluster_contig_df = contig_df_from_cluster_dict(cluster_dict)

            dims = [i for i in cluster_contig_df.columns if 'dim' in i]
            min_dims = 2 * len(dims)
            if min_dims > cluster_contig_df['contig'].size:
                min_dims = 2
            if not pk:
                pk = int(np.log(cluster_contig_df['contig'].size))
                if pk > 15:
                    pk = 15
            if pk < min_dims:
                pk = min_dims

            if cluster_mode == 'HDBSCAN' or not cluster_mode:
                logging.debug(
                    'Trying with HDBSCAN with depth next for cluster {0}.'.format(shorten_cluster_names(cluster)))
                hdbscan_out = hdbscan_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan,
                                                   depth=include_depth, hdbscan_epsilon=hdbscan_epsilon,
                                                   hdbscan_min_samples=hdbscan_min_samples, dist_metric=dist_metric)
                if hdbscan_out[0]:
                    outputs.append(hdbscan_out)
                    continue

                outputs.append([{}, cluster])
                continue
        else:
            if clust_pur == 0:
                logging.debug('Could not calculate purity for cluster {0}.'
                              ' Leaving at 0 and skipping.'.format(shorten_cluster_names(cluster)))
                outputs.append([{}, cluster])
                continue
            elif clust_comp < completeness_threshold:
                logging.debug('Cluster {0} below completeness threshold with'
                              ' {1} and purity {2}. Skipping.'.format(shorten_cluster_names(cluster), clust_comp,
                                                                      clust_pur))
                outputs.append([{}, cluster])
                continue
            else:
                if purity_threshold > clust_pur:
                    logging.warning(shorten_cluster_names(cluster), purity_threshold, clust_pur)
                    logging.warning(clust_dat)
                    raise Exception
                logging.debug('Cluster {0} meets purity threshold of {1} with'
                             ' {2} and has completeness of {3}.'.format(shorten_cluster_names(cluster),
                                                                        purity_threshold,
                                                                        clust_pur, clust_comp))
                outputs.append([cluster_dict, cluster])
                continue
    return outputs


def divide_clusters_by_depth(ds_clstr_dict, threads, marker_sets_graph, tigrfam2pfam_data_dict, min_purity=90,
                             min_completeness=50, pk=None, cluster_mode=None, include_depth=False, max_tries=15,
                             hdbscan_epsilon=0.25, hdbscan_min_samples=2, dist_metric='manhattan'):

    min_purity = min_purity / 100
    min_completeness = min_completeness / 100
    dict_cp = ds_clstr_dict.copy()
    n_tries = 1

    cluster_list = [[i] + [len(dict_cp[i]['contigs'])] for i in list(dict_cp.keys())]
    cluster_list.sort(key=lambda i: i[1], reverse=True)

    threads_needed = min(threads, len(cluster_list))
    logging.debug(f'Threads needed: {threads_needed}')
    chunks_to_process = [[] for i in range(threads_needed)]
    for i in cluster_list:
        chunks_to_process.sort(key=lambda i: sum([len(cluster_dict[list(cluster_dict.keys())[0]]['contigs'])
                                                  for cluster_dict in i]))
        chunks_to_process[0].append({i[0]: dict_cp[i[0]]})

    inner_max_threads = int(threads / len(chunks_to_process))
    if inner_max_threads < 1:
        inner_max_threads = 1

    previous_dict_keys = list(dict_cp.keys())

    no_progress = False
    split_contigs = []

    with parallel_backend("loky", inner_max_num_threads=inner_max_threads):
        while n_tries <= max_tries and chunks_to_process and not no_progress:
            n_tries += 1
            sub_clstr_res = Parallel(n_jobs=threads_needed) \
                (delayed(get_sub_clusters)(cluster_dict_list, inner_max_threads, marker_sets_graph,
                                           tigrfam2pfam_data_dict, purity_threshold=min_purity,
                                           completeness_threshold=min_completeness,
                                           pk=pk, cluster_mode=cluster_mode,
                                           include_depth=include_depth,
                                           hdbscan_epsilon=hdbscan_epsilon,
                                           hdbscan_min_samples=hdbscan_min_samples,
                                           dist_metric=dist_metric)
                 for cluster_dict_list in chunks_to_process)

            for outputs in sub_clstr_res:
                for output in outputs:
                    if output:
                        if output[1]:
                            dict_cp.pop(output[1])
                            split_contigs.append(output[1])
                        dict_cp.update(output[0])
            if sorted(previous_dict_keys) == sorted(list(dict_cp.keys())):
                no_progress = True
                continue
            previous_dict_keys = list(dict_cp.keys())
            cluster_list = []
            for res_list in sub_clstr_res:
                for res in res_list:
                    if res:
                        for cluster in list(res[0].keys()):
                            cluster_list.append([cluster, len(res[0][cluster]['contigs'])])
            cluster_list.sort(key=lambda i: i[1], reverse=True)
            chunks_to_process = [[] for i in range(threads)]
            for i in cluster_list:
                chunks_to_process.sort(key=lambda i: sum(
                    [len(cluster_dict[list(cluster_dict.keys())[0]]['contigs']) for cluster_dict in i]))
                chunks_to_process[0].append({i[0]: dict_cp[i[0]]})
    logging.debug('Clustering iterations until end: {0}'.format(n_tries - 1))
    start = timer()
    dict_cp = {k: v for k, v in dict_cp.items() if v.get('purity', 0) >= min_purity
               and v.get('completeness', 0) >= min_completeness}

    end = timer()
    logging.debug('Added purity and completeness stats to the final dict in {0}s.'.format(int(end - start)))
    return dict_cp, split_contigs


def contig_df_from_cluster_dict(cluster_dict):
    clust_contig_df_rows = []
    dims = [i for i in list(cluster_dict[list(cluster_dict.keys())[0]].keys()) if 'dim' in i]
    depths = [i for i in list(cluster_dict[list(cluster_dict.keys())[0]].keys()) if 'depth' in i]
    clust_contig_df_cols = ['contig', 'essential'] + dims + depths + ['cluster', 'purity', 'completeness']
    # clust_contig_df_col_dtypes = ['string[pyarrow]', 'string[pyarrow]'] + ['float32' for i in dims] \
    #                              + ['float32' for i in depths] + ['string[pyarrow]', 'float32', 'float32']
    clust_contig_df_ind = []
    clust_contig_df_ind_init = 0
    for cluster in cluster_dict:
        warning_off = False
        for index, contig in enumerate(cluster_dict[cluster]['contigs']):
            try:
                new_row = [contig, cluster_dict[cluster]['essential'][index]] \
                          + [cluster_dict[cluster][dim][index] for dim in dims] \
                          + [cluster_dict[cluster][depth][index] for depth in depths] \
                          + [cluster, cluster_dict[cluster].get('purity', 0), cluster_dict[cluster].get('completeness', 0)]
            except KeyError:
                if not warning_off:
                    logging.warning('Something went wrong while fetching cluster data for cluster {0}:'
                                    ' {1}'.format(cluster, [cluster_dict[cluster].get('purity'),
                                                            cluster_dict[cluster].get('completeness')]))
                warning_off = True
                new_row = [contig, cluster_dict[cluster]['essential'][index]] \
                          + [cluster_dict[cluster][dim][index] for dim in dims] \
                          + [cluster_dict[cluster][depth][index] for depth in depths] \
                          + [cluster, cluster_dict[cluster].get('purity', 0),
                             cluster_dict[cluster].get('completeness', 0)]
            clust_contig_df_rows.append(new_row)
            clust_contig_df_ind.append(clust_contig_df_ind_init)
            clust_contig_df_ind_init += 1
    clust_contig_df = pd.DataFrame(clust_contig_df_rows, clust_contig_df_ind, clust_contig_df_cols)
    downcast_pd_df(clust_contig_df)
    return clust_contig_df


def write_scatterplot(df, hue, file_path=None):
    palette = sns.color_palette("husl", len(set(hue)))
    sorted_set = []
    for i in hue:
        if i not in sorted_set:
            sorted_set.append(i)
    for index, i in enumerate(sorted_set):
        if i == 'N':
            palette[index] = ('silver')

    dims = [i for i in df.keys() if 'dim' in i]
    if len(dims) > 3 or len(dims) == 1:
        logging.warning('More than 3 or 1 dimensions. Cant plot.')
        return
    elif len(dims) == 3:
        df['cluster'] = hue
        palette_dict = {cluster: color for cluster, color in zip(sorted_set, palette)}
        fig = plt.figure()
        ax = Axes3D(fig)
        x = df['dim1']
        y = df['dim2']
        z = df['dim3']
        ax.scatter(x, y, z, s=1.5, c=df['cluster'].map(palette_dict))
        fake_list = []
        for i, e in zip(sorted_set, palette):
            fake_dot = ax.scatter(df['dim1'].iloc[0], df['dim2'].iloc[0], df['dim3'].iloc[0], s=1.5, color=e)
            fake_list.append(fake_dot)
        ax.legend(fake_list, [str(i) for i in sorted_set], markerscale=2)
        plt.show()
        if file_path:
            plt.savefig(file_path)
            plt.clf()
    elif len(dims) == 2:
        scatter_plot = sns.scatterplot(data=df, x="dim1", y="dim2", hue=hue, sizes=2.5, s=2., palette=palette)
        scatter_plot.legend(fontsize=3, title="Clusters", title_fontsize=4, ncol=1,
                            bbox_to_anchor=(1.01, 1), borderaxespad=0)
        plt.show()
        if file_path:
            scatter_plot.get_figure().savefig(file_path)
            scatter_plot.get_figure().clf()  # Clear figure


def load_fasta(fasta):
    sequence_dict = {}
    with open(fasta) as f:
        current_header = None
        for line in f:
            # Get header
            if line.startswith('>'):
                # Transform sublists of sequences from previous header to one list.
                if current_header:
                    sequence_dict[current_header] = ''.join(sequence_dict[current_header])
                line = line.strip('\n >')
                sequence_dict[line] = []
                current_header = line
            # Get sequences
            else:
                line = line.strip('\n ')
                sequence_dict[current_header].append(line)
        # Transform sublists of sequences from last header to one list.
        sequence_dict[current_header] = ''.join(sequence_dict[current_header])
    return sequence_dict


def shorten_cluster_names(cluster):
    new_name = []
    sub_clusters = cluster.split('.')
    previous_sub_clstr = None
    for index, sub_clstr in enumerate(sub_clusters):
        if sub_clstr == previous_sub_clstr:
            continue
        previous_sub_clstr = sub_clstr
        counter = 0
        while index + counter + 1 <= len(sub_clusters) - 1 and sub_clstr == sub_clusters[index + counter + 1]:
            counter += 1
        if counter > 0:
            new_name.append(sub_clstr + 'x' + str(counter + 1))
        else:
            new_name.append(sub_clstr)
    new_name = '.'.join(new_name)
    return new_name


def write_bins(cluster_dict, assembly, min_comp=40, min_pur=90, bin_dir='bins'):
    assembly_dict = load_fasta(assembly)
    # Create bin folder, if it doesnt exist.
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    for cluster in cluster_dict:
        if cluster_dict[cluster]['purity'] >= min_pur / 100 and cluster_dict[cluster]['completeness'] >= min_comp / 100:
            new_cluster_name = shorten_cluster_names(cluster)
            if re.match(r"I[0-9]+R[0-9]+\.[0-9]+", cluster):
                new_cluster_name = re.split(r'\D+', new_cluster_name)[1:]
                new_cluster_name = 'I{0}R{1}.{2}'.format('%02.d' % (int(new_cluster_name[0])),
                                                         '%02.d' % (int(new_cluster_name[1])),
                                                         '.'.join(['%06.d' % (int(e)) for e in new_cluster_name[2:]]))
            bin_name = '_'.join(['binny']
                                + [new_cluster_name]
                                + ['C' + str(int(round(cluster_dict[cluster]['completeness'] * 100, 0)))]
                                + ['P' + str(int(round(cluster_dict[cluster]['purity'] * 100, 0)))]
                                + [cluster_dict[cluster]['taxon'].replace(' ', '_')])
            bin_file_name = bin_name + '.fasta'
            bin_out_path = bin_dir / bin_file_name
            with open(bin_out_path, 'w') as out_file:
                for contig in cluster_dict[cluster]['contigs']:
                    try:
                        out_file.write('>' + contig + '\n' + assembly_dict.get(contig) + '\n')
                    except TypeError:
                        logging.warning(f'{contig} in cluster_dict oif cluster {cluster} not in assembly_dict')
                        logging.warning('{0}'.format(cluster_dict[cluster]['contigs']))


########################################################################################################################

def rec_comp(kmer):
    base_pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc_kmer = "".join(base_pairs.get(base, base) for base in reversed(kmer))
    return rc_kmer


def kmer_is_canonical(kmer):
    rc_kmer = rec_comp(kmer)
    if kmer <= rc_kmer:
        canonical_kmer = kmer
        is_canonical = True
    else:
        canonical_kmer = rc_kmer
        is_canonical = False
    return canonical_kmer, is_canonical


def frequencyCount(string, substr):
    count = 0
    pos = 0
    while (True):
        pos = string.find(substr, pos)
        if pos > -1:
            count = count + 1
            pos += 1
        else:
            break
    return count


def kmer_counter(kmer_list_can, sequence):
    kmer_count_raw = [frequencyCount(sequence, kl_km) + frequencyCount(sequence, rec_comp(kl_km))
                      if kl_km != rec_comp(kl_km) else frequencyCount(sequence, kl_km) for kl_km in kmer_list_can]
    sum_can_kmers = sum(kmer_count_raw)
    kmer_freq = [kmer_count / sum_can_kmers if kmer_count > 0 else 0 for kmer_count in kmer_count_raw]
    return kmer_freq


def get_contig_kmer_freq(can_k_mers, assembly_chunk):
    chunk_list = []
    for contig in assembly_chunk:
        start = timer()
        kmer_count_list = kmer_counter(can_k_mers, contig[1])
        end = timer()
        if end - start >= 30:
            logging.debug('Counting k-mers in {0} took longer than 30s with {1}s.'.format(contig[0], int(end - start)))
        logging.debug('Counted k-mers in {0} in {1}s.'.format(contig[0], end - start))
        contig_kfreq = [contig[0]] + kmer_count_list
        chunk_list.append(contig_kfreq)
    return chunk_list


def get_contig_kmer_matrix(contig_list, ksize_list, n_jobs=1):
    contig_kmer_freq_matrix = []

    contig_list = [i + [len(i[1])] for i in contig_list]
    contig_list.sort(key=lambda i: i[2], reverse=True)
    start = timer()
    chunks_to_process = [[] for i in range(n_jobs)]

    list_pos = 0
    for contig in contig_list:
        chunks_to_process[list_pos].append(contig)
        list_pos += 1
        if list_pos + 1 > len(chunks_to_process):
            list_pos = 0
    end = timer()
    logging.debug('Created load balanced list in {0}s.'.format(int(end - start)))
    # Try to free mem
    del contig_list
    n_rounds = 0
    logging.info('Calculating k-mer frequencies of sizes: {0}.'.format(', '.join([str(i) for i in ksize_list])))
    for ksize in ksize_list:
        d = {i: ['A', 'T', 'G', 'C'] for i in range(ksize)}
        kmer_list = [''.join(combination) for combination in itertools.product(*[d[k] for k in sorted(d.keys())])]
        kmer_list_can = sorted(list({i for i in kmer_list if kmer_is_canonical(i)[1]}))
        if n_rounds == 0:
            contig_kmer_freq_matrix.append(['contig'] + kmer_list_can)
        else:
            contig_kmer_freq_matrix[0] = contig_kmer_freq_matrix[0] + kmer_list_can
        with parallel_backend("loky"):
            contig_kmer_freq_matrix_chunks = Parallel(n_jobs=n_jobs) \
                (delayed(get_contig_kmer_freq)(kmer_list_can, chunks)
                 for chunks in chunks_to_process)
        contig_counter = 0
        for chunk in contig_kmer_freq_matrix_chunks:
            for contig_freq in chunk:
                contig_counter += 1
                if n_rounds == 0:
                    contig_kmer_freq_matrix.append(contig_freq)
                else:
                    if contig_kmer_freq_matrix[contig_counter][0] != contig_freq[0]:
                        logging.warning(contig_kmer_freq_matrix[contig_counter][0], contig_freq[0])
                        raise Exception
                    contig_kmer_freq_matrix[contig_counter] = contig_kmer_freq_matrix[contig_counter] + contig_freq[1:]
        logging.debug('Finished counting k-mer frequencies for size {0}.'.format(ksize))
        n_rounds += 1
    return contig_kmer_freq_matrix


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def get_perp(n_data_points):
    return int(np.log(n_data_points))


def gff2low_comp_feature_dict(annotation_file):
    contig_dict = {}
    # Get all info for each line
    with open(annotation_file, 'r') as af:
        for line in af:
            line = line.strip(' \t\n').split('\t')
            contig, feature_start, feature_end, attributes = line[0], line[3], line[4], line[8].split(';')
            # Check if feature has no start and/or stop position or length = 0
            if not feature_start or not feature_end or int(feature_start) - int(feature_end) == 0:
                continue
            # Only add line if feature was annotated as essential gene
            for attribute in attributes:
                attribute_name, attribute_id = attribute.split('=')[0], attribute.split('=')[1]
                if attribute_name == 'product' and 'ribosomal RNA' in attribute_id:
                    if not contig_dict.get(contig):
                        contig_dict[contig] = [(int(feature_start), int(feature_end))]
                    else:
                        contig_dict[contig].append((int(feature_start), int(feature_end)))
                elif attribute_name == 'rpt_family':
                    if not contig_dict.get(contig):
                        contig_dict[contig] = [(int(feature_start), int(feature_end))]
                    else:
                        contig_dict[contig].append((int(feature_start), int(feature_end)))
    return contig_dict


def mask_rep_featrues(contig_rep_feature_dict, contig_list):
    for index, contig in enumerate(contig_list):
        if contig_rep_feature_dict.get(contig[0]):
            regions = contig_rep_feature_dict.get(contig[0])
            for region in regions:
                if region[0] - 1 <= 0:
                    start_region = 0
                    start2feature = 0
                else:
                    start_region = region[0] - 1
                    start2feature = region[0] - 2
                contig_list[index][1] = contig[1][:start2feature] + contig[1][start_region:(region[1] - 1)].lower() \
                                        + contig[1][region[1]:]


def downcast_pd_df(df):
    for column in df:
        if df[column].dtype == 'float64':
            df[column] = pd.to_numeric(df[column], downcast='float')
        if df[column].dtype == 'int64':
            df[column] = pd.to_numeric(df[column], downcast='integer')
        if df[column].dtype == 'object':
            df[column] = df[column].astype('string[pyarrow]')


def get_initial_clusters(contig_df, cluster_alg, threads, include_depth, hdbscan_epsilon, hdbscan_min_samples,
                         dist_metric, marker_sets_graph, tigrfam2pfam_data_dict):

    print("get_initial_clusters, contig_df['contig']:", len(contig_df['contig'].to_list()))
    print("get_initial_clusters, contig_df['contig']:", len(set(contig_df['contig'].to_list())))

    if len(contig_df['contig'].to_list()) != len(set(contig_df['contig'].to_list())):
        raise Exception

    init_hdbscan_epsilon = hdbscan_epsilon
    init_hdbscan_min_samples = hdbscan_min_samples
    init_pk = get_perp(contig_df['contig'].size)
    dims = [i for i in contig_df.columns if 'dim' in i and '_' not in i]
    if init_pk < len(dims) * 2:
        init_pk = len(dims) * 2
    init_clust_dict = {}
    # overfull_clusters = True
    ic_attempts = 0
    cluster_counter = 0
    max_init_clstr_tries = 5
    working_df = contig_df
    logging.info('Running initial clustering.')
    while not working_df.empty and ic_attempts <= max_init_clstr_tries:
        if ic_attempts > 0:
            logging.debug('Splitting highly contaminated clusters again with adjusted parameters')
        # Get some clusters
        clust_dict, labels = run_initial_scan(working_df, cluster_alg, threads, include_depth=include_depth, pk=init_pk,
                                              hdbscan_epsilon=init_hdbscan_epsilon, hdbscan_min_samples=init_hdbscan_min_samples,
                                              dist_metric=dist_metric)

        if ic_attempts > 0:
            logging.debug('New clusters after adjusting: {}'.format(len(clust_dict.keys())))
        if ic_attempts == 0:
            initial_clusters = len(clust_dict.keys())
        working_df = pd.DataFrame()
        # Check if clusters are massively contaminated and thus might take a long time to assess
        # with divide_clusters_by_depth
        for cluster in list(clust_dict.keys()):
            n_markers = len([marker for marker_list in clust_dict[cluster]['essential'] for marker in marker_list.split(',')])
            if len(clust_dict.keys()) < 10:
                logging.debug('cluster: {}, n_contigs: {}, n_markers: {}.'.format(cluster, len(clust_dict[cluster]['contigs']),
                                                                                 n_markers))
            if n_markers > 20000:
                cluster_essential_genes = [gene for genes in clust_dict.get(cluster, {}).get('essential')
                                           for gene in genes.split(',') if not gene == 'non_essential']
                marker_set = choose_checkm_marker_set(cluster_essential_genes, marker_sets_graph,
                                                      tigrfam2pfam_data_dict, max_depth_lvl=0)
                taxon, comp, pur = marker_set[0], marker_set[1], marker_set[2]

                if comp >= 0.80 and pur <= 0.20:
                    if comp == 1.00:
                        comp = '>= 1.00'
                    logging.debug(f'Cluster {cluster} is massive and heavily contaminated with completeness {comp} '
                                  f'and purity {pur}.')
                    # Set contaminated large clusters as working_df
                    working_df = contig_df[contig_df['contig'].isin(clust_dict[cluster]['contigs'])]
                else:
                    new_clust_id = str(cluster_counter)
                    init_clust_dict[new_clust_id] = clust_dict.pop(cluster)
                    cluster_counter += 1
            else:
                new_clust_id = str(cluster_counter)
                init_clust_dict[new_clust_id] = clust_dict.pop(cluster)
                cluster_counter += 1
        # Adjust clustering params
        init_hdbscan_epsilon = max(0, init_hdbscan_epsilon - 0.125)
        init_hdbscan_min_samples = hdbscan_min_samples - 1 if hdbscan_min_samples > 1 else None
        init_pk = max(2, init_pk - 2)
        cluster_counter += 1
        ic_attempts += 1
    if ic_attempts > 1:
        # Add leftover unsplittable clusters back
        if not working_df.empty and ic_attempts > max_init_clstr_tries:
            logging.debug(f'Could not split contaminated clusters further with {max_init_clstr_tries} attempts.'
                         f' Passing them on.')
            for cluster in list(clust_dict.keys()):
                new_clust_id = str(cluster_counter)
                init_clust_dict[new_clust_id] = clust_dict.pop(cluster)
                cluster_counter += 1
        logging.debug(f'Initially created {initial_clusters} clusters.')
    logging.info('Final number of clusters: {}.'.format(len(init_clust_dict.keys())))

    return init_clust_dict


def binny_iterate(contig_data_df, threads, marker_sets_graph, tigrfam2pfam_data_dict, min_purity, min_completeness,
                  max_iterations=4, embedding_iteration=1, max_tries=2, include_depth_initial=False,
                  include_depth_main=True, hdbscan_epsilon=0.25, hdbscan_min_samples_range=(2, 3, 4, 5),
                  dist_metric='manhattan', contigs2clusters_out_path='intermediary'):
    leftovers_df = contig_data_df.copy()

    n_iterations = 1
    good_clusters = {}
    all_binned = False

    hdbs_min_samp_ind = 0
    # if embedding_iteration == 1:
    #     logging.info('Running first round with minimum purity of 95% and extended number of iterations.')
    #     hdbscan_min_samples_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50]
    #     max_iterations = len(hdbscan_min_samples_range)
    #     min_purity = 95

    # if embedding_iteration % 2 != 0 and embedding_iteration > 1:
    #     hdbscan_min_samples_range = hdbscan_min_samples_range[::-1]

    good_bins_per_round_list = []

    while n_iterations <= max_iterations:
        logging.info(f'Clustering round {n_iterations}.')
        hdbscan_min_samples = hdbscan_min_samples_range[hdbs_min_samp_ind]
        hdbs_min_samp_ind += 1
        if hdbs_min_samp_ind == len(hdbscan_min_samples_range):
            hdbs_min_samp_ind = 0

        init_clust_dict = get_initial_clusters(leftovers_df, 'HDBSCAN', threads, include_depth_initial,
                                                       hdbscan_epsilon, hdbscan_min_samples, dist_metric,
                                                       marker_sets_graph, tigrfam2pfam_data_dict)

        if n_iterations == 1:
            for cluster in init_clust_dict:
                init_clust_dict[cluster]['purity'] = 0
                init_clust_dict[cluster]['completeness'] = 0
                init_clust_dict[cluster]['taxon'] = 'none'

        logging.debug('Attempting to find sub-clusters using HDBSCAN.')
        new_clust_dict, split_clusters = divide_clusters_by_depth(init_clust_dict, threads, marker_sets_graph,
                                                                  tigrfam2pfam_data_dict, int(min_purity),
                                                                  int(min_completeness), cluster_mode='HDBSCAN',
                                                                  include_depth=include_depth_main, max_tries=max_tries,
                                                                  hdbscan_epsilon=hdbscan_epsilon,
                                                                  hdbscan_min_samples=hdbscan_min_samples,
                                                                  dist_metric=dist_metric)


        new_clust_dict = {'I{0}R{1}.'.format(embedding_iteration, n_iterations) + k: v for k, v in
                          new_clust_dict.items()}
        n_good_bins_round = len(new_clust_dict.keys())
        good_bins_per_round_list.append(n_good_bins_round)
        logging.debug(f'Good bins this round: {n_good_bins_round}.')
        iteration_clust_dict = {**new_clust_dict, **init_clust_dict}

        good_clusters.update(new_clust_dict)
        if list(good_clusters.keys()):
            good_clust_contig_df = contig_df_from_cluster_dict(good_clusters)

        iteration_clust_df = cluster_df_from_dict(iteration_clust_dict)
        iteration_clust_df.to_csv('{0}/iteration_{1}_contigs2clusters.tsv'.format(contigs2clusters_out_path,
                                                                                  n_iterations), sep='\t', index=False)
        iteration_clust_contig_df = contig_df_from_cluster_dict(iteration_clust_dict)
        downcast_pd_df(iteration_clust_contig_df)

        # Plot sub-clustering
        conditions = [(iteration_clust_contig_df['purity'] >= min_purity / 100) & (
                iteration_clust_contig_df['completeness'] >= min_completeness / 100),
                      (iteration_clust_contig_df['purity'] < min_purity / 100) | (
                              iteration_clust_contig_df['completeness'] < min_completeness / 100)]
        values = [iteration_clust_contig_df['cluster'], 'N']
        iteration_clust_contig_df['above_thresh'] = np.select(conditions, values)
        iteration_clust_contig_df = contig_data_df.merge(iteration_clust_contig_df, how='outer', on='contig',
                                                         suffixes=(None, '_y'))
        iteration_clust_contig_df = iteration_clust_contig_df[
            iteration_clust_contig_df.columns.drop(list(iteration_clust_contig_df.filter(regex='_y')))]
        iteration_clust_contig_df['above_thresh'] = iteration_clust_contig_df['above_thresh'].fillna('N')

        iteration_clust_contig_df = iteration_clust_contig_df[
            ~((iteration_clust_contig_df['contig'].duplicated(keep=False))
              & (iteration_clust_contig_df['above_thresh'] == 'N'))]

        leftovers_df = iteration_clust_contig_df[iteration_clust_contig_df['above_thresh'] == 'N']
        if list(good_clusters.keys()):
            leftovers_df = leftovers_df[~leftovers_df['contig'].isin(good_clust_contig_df['contig'].tolist())]
        n_iterations += 1
        if len(leftovers_df.index) == 0:
            logging.info('No leftover contigs to cluster. Exiting.')
            all_binned = True
            break
    downcast_pd_df(iteration_clust_contig_df)

    good_bins_per_round_list.sort()
    mid = len(good_bins_per_round_list) // 2
    median_bins_per_round = (good_bins_per_round_list[mid] + good_bins_per_round_list[~mid]) / 2

    # Cheat a bit on first round
    if embedding_iteration == 1:
        if sum(good_bins_per_round_list) >= 4:
            median_bins_per_round = 2

    return good_clusters, init_clust_dict, all_binned, median_bins_per_round


def get_single_contig_bins(essential_gene_df, good_bins_dict, n_dims, marker_sets_graph, tigrfam2pfam_data_dict,
                           threads=1):
    essential_gene_lol = essential_gene_df.values.tolist()
    cluster_list = [[i[0], i[1].split(',')] for i in essential_gene_lol if len(set(i[1].split(','))) >= 40]
    cluster_list.sort(key=lambda i: i[1], reverse=True)
    start = timer()
    chunks_to_process = [[] for i in range(threads)]
    for i in cluster_list:
        chunks_to_process.sort(
            key=lambda i: sum([len(cont_ess_gene_data[1]) for cont_ess_gene_data in i]))
        chunks_to_process[0].append(i)

    with parallel_backend("loky"):
        single_contig_bin_dict_list = Parallel(n_jobs=threads) \
            (delayed(asses_contig_completeness_purity)
             (contig_data, n_dims, marker_sets_graph,
              tigrfam2pfam_data_dict)
             for contig_data in chunks_to_process)

    for bin_dict_sub_list in single_contig_bin_dict_list:
        for bin_dict in bin_dict_sub_list:
            good_bins_dict.update(bin_dict)
    end = timer()
    logging.debug('Finished searching for single contig bins in {0}s.'.format(int(end - start)))
    return list(good_bins_dict.keys())


def asses_contig_completeness_purity(essential_gene_lol, n_dims, marker_sets_graph, tigrfam2pfam_data_dict):
    single_contig_bins = []
    for contig_data in essential_gene_lol:
        all_ess = contig_data[1]
        marker_set = choose_checkm_marker_set(all_ess, marker_sets_graph, tigrfam2pfam_data_dict, max_depth_lvl=4)
        taxon, comp, pur = marker_set[0], marker_set[1], marker_set[2]
        if pur > 0.90 and comp > 0.925:
            bin_dict = {contig_data[0]: {'depth1': np.array([None]), 'contigs': np.array([contig_data[0]]),
                                         'essential': np.array(all_ess), 'purity': pur, 'completeness': comp,
                                         'taxon': taxon}}
            for dim in range(n_dims):
                bin_dict[contig_data[0]]['dim' + str(dim + 1)] = np.array(np.array([None]))
            single_contig_bins.append(bin_dict)
    return single_contig_bins


def tigrfam2pfam_dict(tigrfam2pfam_file):
    tigrfam2pfam_dict = {}
    with open(tigrfam2pfam_file, 'r') as f:
        for line in f:
            line = line.strip('\n \t').split('\t')
            if not tigrfam2pfam_dict.get(line[1]):
                tigrfam2pfam_dict[line[1]] = [line[0].replace('pfam', 'PF')]
            else:
                tigrfam2pfam_dict[line[1]].append(line[0].replace('pfam', 'PF'))
            if not tigrfam2pfam_dict.get(line[0]):
                tigrfam2pfam_dict[line[0].replace('pfam', 'PF')] = [line[1]]
            else:
                tigrfam2pfam_dict[line[0].replace('pfam', 'PF')].append(line[1])
    return tigrfam2pfam_dict


def load_checkm_markers(marker_file):
    linegae_dict = {0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species'}
    tms_data = nx.DiGraph()
    with open(marker_file, 'r') as f:
        # next(f)
        for line in f:
            line = line.strip('\n \t').split('\t')
            if (int(line[5]) < 58
                    or 'Tenericutes' in line[2]
                    or 'Haemophilus' in line[2]
                    # or 'Enterobacteriales' in line[2]
                    or 'Rickettsiales' in line[2]
                    or 'Chlamydiae' in line[2]):
                logging.debug('Skipping loading of markers from {0}.'.format(line[1]))
                continue
            # marker_sets = line[-1].replace(']), set([', ';')
            # for char in '[](\'() \"':
            #     marker_sets = marker_sets.replace(char, '').replace('set', '')
            # marker_sets = [[marker.split('.')[0] for marker in marker_set.split(',')] for marker_set in
            #                marker_sets.split(';') if not marker_set.split(',') == ['']]
            marker_sets = eval(line[-1])
            marker_sets = [{marker.split('.')[0] for marker in marker_set} for marker_set in
                           marker_sets]
            lineage = line[2].split(';')
            lineage = [tax + '[' + linegae_dict[tax_ind] + ']' if tax == lineage[tax_ind - 1] and len(lineage) > 1
                       else tax for tax_ind, tax in enumerate(lineage)]
            if len(lineage) == 7:
                lineage[-1] = ' '.join(lineage[-2:])

            if len(marker_sets) >= 10:
                tms_data.add_nodes_from(lineage)
                tms_data.add_edges_from([(i, lineage[index + 1]) for index, i in enumerate(lineage)
                                         if not index + 1 == len(lineage)])
                tms_data.nodes[lineage[-1]]['markers'] = int(line[4])
                tms_data.nodes[lineage[-1]]['marker_groups'] = int(line[5])
                tms_data.nodes[lineage[-1]]['marker_sets'] = marker_sets

                # Check if higher level tax marker sets are identical
                if not len(lineage) == 1:
                    for level in lineage[:-1]:
                        try:
                            if tms_data.nodes[level]['marker_sets'] == marker_sets:
                                tms_data.nodes[lineage[-1]]['marker_sets'] = 'is_' + level
                                break
                        except KeyError:
                            logging.warning(lineage)
                            logging.warning(tms_data.nodes[lineage[-1]])
                            raise KeyError
    return tms_data


def get_marker_set_quality(marker_set, marker_list, tigrfam2pfam_data_dict):
    marker_set_t2p_markers_found = [t2p_marker for marker in marker_list
                         for t2p_marker in [marker] + tigrfam2pfam_data_dict.get(marker, [])
                         if t2p_marker in marker_set]

    if not marker_set_t2p_markers_found:
        return [0, '']

    node_markers_list_t2p = []
    for marker in set(marker_set_t2p_markers_found):
        node_markers_list_t2p += tigrfam2pfam_data_dict.get(marker, [])

    node_t2p_markers_list_set = set(marker_set_t2p_markers_found).union(node_markers_list_t2p)

    marker_set_markers_found = [marker for marker in marker_list
                                if marker in marker_set]
    node_markers_list_set = set(marker_set_markers_found)

    if not marker_set_markers_found:
        marker_set_markers_found = marker_set_t2p_markers_found
        node_markers_list_set = set(node_t2p_markers_list_set).intersection(set(marker_set_t2p_markers_found))

    n_set_markers_found = marker_set.intersection(node_t2p_markers_list_set)

    marker_set_completeness = round(len(n_set_markers_found) / len(marker_set), 3)

    if marker_set_completeness > 1:
        marker_set_completeness = 1

    marker_set_marker_purities = [round(1 / marker_set_markers_found.count(marker), 3)
                                  # for marker in set(node_t2p_markers_list_set).intersection(set(marker_set_t2p_markers_found))]
                                  for marker in node_markers_list_set]
    try:
        marker_set_average_purity = round(sum(marker_set_marker_purities)
                                          / len(marker_set_marker_purities), 3)
    except ZeroDivisionError:
        print(marker_set_markers_found)
        print(node_markers_list_set)
        print(marker_set_marker_purities)
        print(marker_set_t2p_markers_found)
        raise ZeroDivisionError

    return [marker_set_completeness, marker_set_average_purity]


def get_marker_list_node_quality(marker_list, node, marker_sets_graph, tigrfam2pfam_data_dict):
    # node_marker_sets = [set(marker_set) for marker_set in marker_sets_graph.nodes.data()[node]['marker_sets']]
    node_marker_sets = marker_sets_graph.nodes.data()[node]['marker_sets']

    node_marker_sets_completenesses = []
    node_marker_sets_purities = []

    for marker_set in node_marker_sets:
        marker_set_stats = get_marker_set_quality(marker_set, marker_list, tigrfam2pfam_data_dict)
        if marker_set_stats:
            node_marker_sets_completenesses.append(marker_set_stats[0])
            if marker_set_stats[1]:
                node_marker_sets_purities.append(marker_set_stats[1])
        else:
            node_marker_sets_completenesses.append(0)

    node_marker_set_completeness = round(sum(node_marker_sets_completenesses)
                                         / len(node_marker_sets_completenesses), 3)
    if node_marker_sets_purities:
        node_marker_set_purity = round(sum(node_marker_sets_purities)
                                       / len(node_marker_sets_purities), 3)
    else:
        node_marker_set_purity = 0

    if node_marker_set_completeness > 1:
        logging.error('Completeness of for marker set {0} is > 1 with {1} for'
                      ' marker list {2}'.format(node, node_marker_set_completeness,
                                                marker_list))
        raise Exception

    return [node_marker_set_completeness, node_marker_set_purity]


def compare_marker_set_stats(marker_set, current_best_marker_set, completeness_variability):
    if marker_set[1] >= current_best_marker_set[1] * completeness_variability:
        current_best_marker_set = marker_set
    return current_best_marker_set


def choose_checkm_marker_set(marker_list, marker_sets_graph, tigrfam2pfam_data_dict, max_depth_lvl=6):
    nodes = [n for n, d in marker_sets_graph.in_degree() if d == 0]
    current_node = nodes[0]
    previous_nodes = None
    best_marker_set = []
    depth_grace_count = 0
    king_lvl_stats = None
    # 0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species'
    current_depth_level = 0
    while list(marker_sets_graph[current_node]) and depth_grace_count < 2 and current_depth_level <= max_depth_lvl:
        current_level_best_marker_set = []
        if previous_nodes == nodes:
            depth_grace_count += 1
            nodes = [sub_node for node in nodes for sub_node in list(marker_sets_graph[node])]
        previous_nodes = nodes
        for index, node in enumerate(nodes):
            if isinstance(marker_sets_graph.nodes.data()[node]['marker_sets'], str):
                logging.debug('Marker set of {0} identical to higher level set {1}.'
                              ' Skipping.'.format(node,
                                                  marker_sets_graph.nodes.data()[node]['marker_sets'].split('_')[
                                                      1]))
                if not best_marker_set:
                    best_marker_set = [node, 0.00, 0.00, 0.00, current_depth_level]
                if not current_level_best_marker_set:
                    current_level_best_marker_set = [node, 0.00, 0.00, 0.00, current_depth_level]
                continue

            node_n_markers = marker_sets_graph.nodes.data()[node]['markers']
            node_n_marker_sets = marker_sets_graph.nodes.data()[node]['marker_groups']
            node_stats = get_marker_list_node_quality(marker_list, node, marker_sets_graph,
                                                      tigrfam2pfam_data_dict)
            if not node_stats:
                continue

            node_marker_set_completeness = node_stats[0]
            node_marker_set_purity = node_stats[1]

            # # Minimize use of Kingdom level sets, which are most error prone
            # if node in ['Bacteria', 'Archaea']:
                # node_marker_set_completeness -= 0.05
                # node_marker_set_purity -= 0.05

            node_marker_set_completeness_score = round(node_marker_set_completeness
                                                       * node_n_marker_sets / node_n_markers * 100, 3)

            if king_lvl_stats: # and node not in ['Bacteria', 'Archaea']:
                current_marker_set = [node, ((node_marker_set_completeness + king_lvl_stats[0]) / 2),
                                      ((node_marker_set_purity + king_lvl_stats[1]) / 2),
                                      node_marker_set_completeness_score, current_depth_level]
            else:
                current_marker_set = [node, node_marker_set_completeness, node_marker_set_purity,
                                      node_marker_set_completeness_score, current_depth_level]

            if not best_marker_set: # or (best_marker_set[0] in ['Bacteria', 'Archaea'] and current_depth_level > 0):
                best_marker_set = [node, node_marker_set_completeness, node_marker_set_purity,
                                   node_marker_set_completeness_score, current_depth_level]
            else:
                # Check if comparing same level sets, if so dont give completeness leeway
                if current_marker_set[-1] == best_marker_set[-1]:
                    completeness_variability = 1.0
                else:
                    completeness_variability = 0.975
                best_marker_set = compare_marker_set_stats(current_marker_set,
                                                              best_marker_set, completeness_variability)

            if not current_level_best_marker_set:  #  or (current_level_best_marker_set[0] in ['Bacteria', 'Archaea']
                                                   # and current_depth_level > 0):
                current_level_best_marker_set = [node, node_marker_set_completeness, node_marker_set_purity,
                                                 node_marker_set_completeness_score, current_depth_level]
            else:
                if current_marker_set[-1] == current_level_best_marker_set[-1]:
                    completeness_variability = 1.0
                else:
                    completeness_variability = 0.975
                current_level_best_marker_set = compare_marker_set_stats(current_marker_set,
                                                                         current_level_best_marker_set,
                                                                         completeness_variability)

        nodes = list(marker_sets_graph[current_level_best_marker_set[0]])

        current_node = current_level_best_marker_set[0]

        if current_level_best_marker_set[0] in ['Bacteria', 'Archaea']:
            king_lvl_stats = [node_marker_set_completeness, node_marker_set_purity]

        current_depth_level += 1

    if best_marker_set:
        return best_marker_set
    else:
        logging.debug('Something went wrong while choosing the best marker set. Markers:'
                      ' {0}; unique: {1}; total {2}.'.format(set(marker_list),
                                                             len(set(marker_list)), len(marker_list)))
        return ['None', 0, 0, 0]


def gold_standard_stats(cont2gen_gs_dict, assembly_dict):
    gs_stats_dict = {}
    for contig in cont2gen_gs_dict:
        cont_genome = cont2gen_gs_dict[contig]
        if not cont2gen_gs_dict[contig] in gs_stats_dict:
            gs_stats_dict[cont_genome] = np.array([len(assembly_dict[contig]), 1])
        else:
            gs_stats_dict[cont_genome] += np.array([len(assembly_dict[contig]), 1])
    return gs_stats_dict


def load_depth_dict(mg_depth_file):
    depth_dict = {}
    with open(mg_depth_file, 'r') as f:
        for line in f:
            line = line.strip('\n \t').split('\t')
            try:
                depth_dict[line[0]] = np.array([float(depth_val) for depth_val in line[1:]])
            except ValueError:
                continue
    return depth_dict


def parse_mantis_cons_annot(mantis_out_annot):
    mantis_data = {}
    with open(mantis_out_annot, 'r') as f:
        next(f)
        for line in f:
            line = line.strip('\n \t').split('\t')
            gene, markers = line[0], [marker.split(':')[-1] for marker in line[6:] if
                                      not marker.split(':')[-1].startswith('DUF')]
            mantis_data[gene] = markers
    return mantis_data


def checkm_hmmer_search2prokka_gff(hmm_checkm_marker_out, prokka_gff, gff_out_path='intermediary'):
    prokka_checkm_gff = '{0}/annotation_CDS_RNA_hmms_checkm.gff'.format(gff_out_path)
    checkm_marker_dict = parse_mantis_cons_annot(hmm_checkm_marker_out)
    with open(prokka_checkm_gff, 'w') as pcg:
        with open(prokka_gff, 'r') as pg:
            for line in pg:
                line = line.strip('\n \t')
                line_annots = line.split('\t')[-1].split(';')
                line_gene = line_annots[0].replace('ID=', '')
                if checkm_marker_dict.get(line_gene):
                    pcg.write(line + ';checkm_marker=' + ','.join(checkm_marker_dict[line_gene]) + '\n')
                else:
                    pcg.write(line + '\n')


def check_sustainable_contig_number(contig_list, min_val, assembly_dict, contig_threshold=5e5):
    n_contigs = len([cont for cont in contig_list if len(assembly_dict[cont]) >= min_val])
    while n_contigs > contig_threshold:
        min_val += 5
        n_contigs = len([cont for cont in contig_list if len(assembly_dict[cont]) >= min_val])
    return min_val


def calc_assembly_nx(assembly_dict, scmags, nx_val):
    mini_dict = {id: id for id in scmags}
    assembly_cont_length_array = np.sort(np.array([len(seq) for contig_id, seq in assembly_dict.items()
                                                   if not mini_dict.get(contig_id)]))[::-1]
    nx_size = int(sum(assembly_cont_length_array) * nx_val / 100)
    size_sum = 0
    i = 0
    nx = 0
    while size_sum <= nx_size:
        length = assembly_cont_length_array[i]
        nx = length
        size_sum += length
        i += 1
    return nx


def iterative_embedding(x_contigs, depth_dict, all_good_bins, starting_completeness, min_purity, min_completeness,
                        threads, n_dim, annot_file, mg_depth_file, single_contig_bins, taxon_marker_sets,
                        tigrfam2pfam_data, main_contig_data_dict, assembly_dict, max_contig_threshold=3.0e5,
                        internal_min_marker_cont_size=0, include_depth_initial=False, max_embedding_tries=50,
                        include_depth_main=True, hdbscan_epsilon_range=(0.250, 0.125), hdbscan_min_samples_range=(0, 2, 4, 6, 8),
                        dist_metric='euclidean', contigs2clusters_out_path='intermediary'):
    np.random.seed(0)
    embedding_tries = 1
    internal_completeness = starting_completeness
    final_try_counter = 0
    pk_factor = 1
    hdbscan_epsilon = hdbscan_epsilon_range[0]
    perp_1 = 10
    perp_2 = 100
    while embedding_tries <= max_embedding_tries:
        if embedding_tries == 1:
            internal_min_marker_cont_size = check_sustainable_contig_number(x_contigs, internal_min_marker_cont_size,
                                                                            assembly_dict, max_contig_threshold)
            round_x = [main_contig_data_dict.get(cont) for cont in x_contigs
                       if len(assembly_dict[cont]) >= internal_min_marker_cont_size]
            round_x_contigs = [cont for cont in x_contigs if len(assembly_dict[cont]) >= internal_min_marker_cont_size]
            round_leftovers_contig_list = [cont for cont in x_contigs
                                           if len(assembly_dict[cont]) < internal_min_marker_cont_size]
        else:
            round_leftovers_contig_list_backup = round_leftovers_contig_list.copy()
            internal_min_marker_cont_size = check_sustainable_contig_number(round_leftovers_contig_list,
                                                                            internal_min_marker_cont_size,
                                                                            assembly_dict, max_contig_threshold)
            round_x = [main_contig_data_dict.get(cont) for cont in round_leftovers_contig_list
                       if len(assembly_dict[cont]) >= internal_min_marker_cont_size]
            round_x_contigs = [cont for cont in round_leftovers_contig_list
                               if len(assembly_dict[cont]) >= internal_min_marker_cont_size]
            round_leftovers_contig_list = [cont for cont in round_leftovers_contig_list
                                           if len(assembly_dict[cont]) < internal_min_marker_cont_size]
            while ((max_contig_threshold < len(round_x_contigs) or len(round_x_contigs) < 5)
                   and internal_min_marker_cont_size > 0):
                if len(round_x_contigs) < 5:
                    logging.info('Less than 5 contigs to bin. Deacreasing min length threshold.')
                    internal_min_marker_cont_size -= 125
                    if internal_min_marker_cont_size < 0:
                        internal_min_marker_cont_size = 0
                else:
                    logging.info('More than {0} contigs to bin. Increasing min length threshold.'.format(
                        max_contig_threshold))
                    internal_min_marker_cont_size += 250
                round_x = [main_contig_data_dict.get(cont) for cont in round_leftovers_contig_list_backup
                           if len(assembly_dict[cont]) >= internal_min_marker_cont_size]
                round_x_contigs = [cont for cont in round_leftovers_contig_list_backup
                                   if len(assembly_dict[cont]) >= internal_min_marker_cont_size]
                round_leftovers_contig_list = [cont for cont in round_leftovers_contig_list_backup
                                               if len(assembly_dict[cont]) < internal_min_marker_cont_size]
            del round_leftovers_contig_list_backup

        logging.info('Running with {0} contigs. Filtered {1} contigs using a min contig size of {2} to stay below'
                     ' {3} contigs'.format(len(round_x_contigs), len(round_leftovers_contig_list),
                                           internal_min_marker_cont_size, max_contig_threshold))

        if len(round_x_contigs) != len(round_x):
            logging.warning('Contig feature data length ({0}) doesnt match contig id'
                            ' list length ({1}). Exiting'.format(len(round_x_contigs), len(round_x)))
            raise Exception

        round_x_depth = np.array([depth_dict[contig] for contig in round_x_contigs])

        # Replace zeroes for clr
        round_x = multiplicative_replacement(round_x)
        # Clr transform
        x_scaled = clr(round_x)
        x_scaled = np.concatenate([round_x_depth, x_scaled], axis=1)
        # x_scaled = np.concatenate([round_x_depth, round_x], axis=1)

        # Manifold learning and dimension reduction.
        logging.info('Running manifold learning and dimensionality-reduction.')
        n_pca_tries = 0

        if len(round_x_contigs) > 30:
            n_comp = 30
        else:
            n_comp = len(round_x_contigs) - 1

        pca = PCA(n_components=n_comp, random_state=0)
        transformer = pca.fit(x_scaled)
        sum_var_exp = sum(pca.explained_variance_ratio_)
        while sum_var_exp <= 0.75 and n_pca_tries <= 100 and len(pca.explained_variance_ratio_) <= 70 \
                and not len(round_x_contigs) <= n_comp:
            pca = PCA(n_components=n_comp, random_state=0)
            transformer = pca.fit(x_scaled)
            sum_var_exp = sum(pca.explained_variance_ratio_)
            n_comp += 2
            n_pca_tries += 1
        var_exp = int(round(sum(pca.explained_variance_ratio_), 3) * 100)
        logging.debug(f'PCA stats: Dimensions: {n_comp}; Amount of variation explained: {var_exp}%.')
        x_pca = transformer.transform(x_scaled)

        # perp = perp_range[tsne_perp_ind]
        # tsne_perp_ind += 1
        # if tsne_perp_ind == len(perp_range):
        #     tsne_perp_ind = 0
        if embedding_tries > 1:
            if perp_1 < 20:
                perp_1 += 2
            else:
                perp_1 = 10
            if perp_2 < 130:
                perp_2 += 5
            else:
                perp_2 = 100
        perplexities = [perp_1, perp_2]

        early_exagg = max(4, min(100, int(len(x_pca) * 0.00025)))


        learning_rate = max(2, int(len(x_pca) / early_exagg))  # learning_rate_factor
        logging.debug(f'tSNE learning rate: {learning_rate}, early_exagg: {early_exagg},'
                     f' perplexies: {perplexities}, pk_factor: {pk_factor}')

        logging.info(f'Initializing embedding.')

        affinities_multiscale_mixture = affinity.Multiscale(x_pca, perplexities=perplexities, metric=dist_metric,
                                                            n_jobs=threads, random_state=0, verbose=0)
        init = initialization.pca(x_pca, random_state=0, verbose=0)
        embedding = TSNEEmbedding(init, affinities_multiscale_mixture, negative_gradient_method="fft", n_jobs=threads,
                                  random_state=0, verbose=0)
        opt_cycle_iter = 250
        ee_iter = 0

        # Early exagg
        embedding.optimize(n_iter=3, exaggeration=early_exagg, inplace=True, learning_rate=learning_rate, momentum=0.5,
                           n_jobs=threads, verbose=0)
        KLD_prev = embedding.kl_divergence
        # KLD_DIFF = 1e10
        ee_iter += opt_cycle_iter
        logging.info('EE iteration {0} - KLD: {1}'.format(ee_iter, round(KLD_prev, 4)))

        while 2000 >= ee_iter <= 750 or KLD_DIFF > KLD * 0.01:  # or KLDRC < maxKLDRC:
            embedding.optimize(n_iter=opt_cycle_iter, exaggeration=early_exagg, inplace=True, learning_rate=learning_rate,
                               momentum=0.5, n_jobs=threads, verbose=0)
            KLD = embedding.kl_divergence
            KLD_DIFF = KLD_prev - KLD
            KLDRC = 100 * KLD_DIFF / KLD_prev
            KLD_prev = KLD

            # if KLD_DIFF >= 0:
            #     embedding_prev = embedding
            # else:
            #     embedding = embedding_prev

            ee_iter += opt_cycle_iter
            logging.info('EE iteration {0} - KLD: {1}, KLDRC: {2}%'.format(ee_iter, round(KLD, 4), round(KLDRC, 2)))

        # Main cycle
        opt_cycle_iter = 250
        main_iter = 0

        learning_rate = max(200, min(64e3, int(len(x_pca) * 0.1)))
        logging.info(f'Main iteration learning rate: {learning_rate}')

        while 10000 >= main_iter <= 1000 or KLD_DIFF > KLD * 0.01:
            embedding.optimize(n_iter=opt_cycle_iter, exaggeration=1, inplace=True, learning_rate=learning_rate,
                               momentum=0.8, n_jobs=threads, verbose=0)
            KLD = embedding.kl_divergence
            KLD_DIFF = KLD_prev - KLD
            KLDRC = 100 * KLD_DIFF / KLD_prev

            # if KLD_DIFF >= 0:
            #     embedding_res_prev = embedding.view(np.ndarray)
            #     KLD_FIN = KLD
            # else:
            #     KLD_FIN = KLD_prev

            KLD_prev = KLD

            main_iter += opt_cycle_iter
            logging.info('Main iteration {0} - KLD: {1}, KLDRC: {2}%'.format(main_iter, round(KLD, 4), round(KLDRC, 2)))

        # if KLD_DIFF >= 0:
        #     embedding_result = embedding.view(np.ndarray)
        # else:
        #     embedding_result = embedding_res_prev
        embedding_result = embedding.view(np.ndarray)

        total_iter = ee_iter + main_iter
        logging.info(f'Finished dimensionality-reduction in {total_iter} iterations.'
                     f' Final KLD: {round(KLD, 4)}')

        # Create coordinate df.
        dim_range = [i + 1 for i in range(n_dim)]
        coord_df = pd.DataFrame(data=embedding_result, index=None, columns=['dim' + str(i) for i in dim_range],
                                dtype='float32')
        coord_df['contig'] = round_x_contigs
        coord_df['contig'] = coord_df['contig'].astype('string[pyarrow]')

        # Reorder
        coord_df = coord_df[['contig'] + ['dim' + str(i) for i in dim_range]]

        # Load data
        contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, coord_df, dims=n_dim,
                                                  coords_from_file=False)

        for i in single_contig_bins:
            if i in coord_df['contig'].tolist():
                logging.warning('Single contig bins found in contig df. Something went wrong, exiting.')
                raise Exception

        # Preserve original contig data
        if embedding_tries == 1:
            contig_data_df_org = contig_data_df.copy()

        if embedding_tries > 1:
            if hdbscan_epsilon > hdbscan_epsilon_range[1]:
                hdbscan_epsilon -= 0.125
            else:
                hdbscan_epsilon = hdbscan_epsilon_range[0]

        # Find bins
        good_bins, final_init_clust_dict, all_binned, median_bins_per_round = binny_iterate(contig_data_df, threads,
            taxon_marker_sets, tigrfam2pfam_data, min_purity, internal_completeness, len(hdbscan_min_samples_range),
            embedding_iteration=embedding_tries, max_tries=2, include_depth_initial=include_depth_initial,
            include_depth_main=include_depth_main, hdbscan_epsilon=hdbscan_epsilon,
            hdbscan_min_samples_range=hdbscan_min_samples_range, dist_metric=dist_metric,
            contigs2clusters_out_path=contigs2clusters_out_path)

        logging.info('Good bins this embedding iteration: {0}.'.format(len(good_bins.keys())))

        if median_bins_per_round < 1 and internal_completeness > min_completeness and final_try_counter == 0 \
           and not all_binned:
            internal_completeness -= 10
            if internal_completeness < 90 and min_purity < 90:
                min_purity = 90
            logging.info(f'Median of good bins per round < 1. Minimum completeness lowered to {internal_completeness}.')
        elif median_bins_per_round < 1 and final_try_counter <= 4 and not all_binned \
                and not internal_min_marker_cont_size > prev_round_internal_min_marker_cont_size:
            if final_try_counter == 0:
                max_contig_threshold *= 1.25
            internal_min_marker_cont_size = 2000 - 500 * final_try_counter
            final_try_counter += 1
            logging.info(f'Running with contigs >= {internal_min_marker_cont_size}bp, minimum completeness {internal_completeness}.')
        elif len(list(good_bins.keys())) < 2 and not all_binned:
            logging.info('Reached min completeness and min contig size. Exiting embedding iteration')
            break

        prev_round_internal_min_marker_cont_size = internal_min_marker_cont_size

        round_clust_dict = {**good_bins, **final_init_clust_dict}

        round_clust_contig_df = contig_df_from_cluster_dict(round_clust_dict)

        # Filter round results and ensure no duplicates
        conditions = [(round_clust_contig_df['purity'] >= min_purity / 100) & (
                round_clust_contig_df['completeness'] >= internal_completeness / 100),
                      (round_clust_contig_df['purity'] < min_purity / 100) | (
                              round_clust_contig_df['completeness'] < internal_completeness / 100)]
        values = [round_clust_contig_df['cluster'], 'N']
        round_clust_contig_df['above_thresh'] = np.select(conditions, values)
        round_clust_contig_df = contig_data_df.merge(round_clust_contig_df, how='outer', on='contig',
                                                     suffixes=(None, '_y'))
        round_clust_contig_df = round_clust_contig_df[
            round_clust_contig_df.columns.drop(list(round_clust_contig_df.filter(regex='_y')))]
        round_clust_contig_df['above_thresh'] = round_clust_contig_df['above_thresh'].fillna('N')

        round_leftovers = round_clust_contig_df[round_clust_contig_df['above_thresh'] == 'N']

        if len(set(round_leftovers['contig'].tolist())) != len(round_leftovers['contig'].tolist()):
            logging.debug('{0} dupilcates in leftovers.'
                          ' Removing them.'.format(len(round_leftovers['contig'].tolist())
                                                   - len(set(round_leftovers['contig'].tolist()))))
            round_leftovers.drop_duplicates(subset=['contig'])

        if list(good_bins.keys()):
            good_bin_contig_df = contig_df_from_cluster_dict(good_bins)
            round_leftovers = round_leftovers[~round_leftovers['contig'].isin(good_bin_contig_df['contig'].tolist())]

        round_leftovers_contig_list = round_leftovers_contig_list + round_leftovers['contig'].tolist()

        all_good_bins.update(good_bins)
        embedding_tries += 1

        all_contigs = []
        for bin in all_good_bins:
            all_contigs.extend(all_good_bins[bin]['contigs'])
        if len(all_contigs) != len(set(all_contigs)):
            logging.warning('WARNING: {0} duplicate contigs in bins found! Exiting.'.format(len(all_contigs)
                                                                                            - len(set(all_contigs))))
            raise Exception
        logging.info('Total number of good bins: {0}.'.format(len(all_good_bins.keys())))
        if len(round_leftovers.index) == 0:
            break
    return all_good_bins, contig_data_df_org, min_purity
