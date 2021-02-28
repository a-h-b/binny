"""
Created on Wed Feb 22 10:50:35 2020

@author: oskar.hickl
"""

from pathlib import Path
from joblib import parallel_backend, Parallel, delayed
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from unidip import UniDip


def unify_multi_model_genes(gene):
    # Dict to unify genes with multiple models
    hmm_dict = {'TIGR00388': 'glyS', 'TIGR00389': 'glyS', 'TIGR00471': 'pheT', 'TIGR00472': 'pheT', 'TIGR00408': 'proS',
                'TIGR00409': 'proS', 'TIGR02386': 'rpoC', 'TIGR02387': 'rpoC'}
    # Replace gene with unified name if it has multiple models
    if gene in hmm_dict:
        gene = hmm_dict.get(gene)
    return gene


def gff2ess_gene_df(annotation_file):
    contig_dict = {}
    # Get all info for each line
    with open(annotation_file, 'r') as af:
        for line in af:
            line = line.strip(' \t\n').split('\t')
            contig, contig_start, contig_end, attributes = line[0], line[3], line[4], line[8].split(';')
            # Check if feature has no start and/or stop position or length = 0
            if not contig_start or not contig_end or int(contig_start)-int(contig_end) == 0:
                continue
            # Only add line if feature was annotated as essential gene
            for attribute in attributes:
                attribute_name, attribute_id = attribute.split('=')[0], attribute.split('=')[1]
                if attribute_name == 'essential' and attribute_id:
                    if not contig_dict.get(contig):
                        contig_dict[contig] = [unify_multi_model_genes(attribute_id)]
                    else:
                        contig_dict[contig].append(unify_multi_model_genes(attribute_id))
    # Create list of lists to join dict items as string and build data frame
    annot_cont_line_list = [[contig, ','.join(contig_dict.get(contig))] for contig in contig_dict]
    annot_cont_df = pd.DataFrame(annot_cont_line_list, columns=['contig', 'essential'])
    return annot_cont_df


def load_and_merge_cont_data(annot_file, depth_file, vizbin_coords):
    print('Reading Prokka annotation.')
    annot_df = gff2ess_gene_df(annot_file)
    print('Reading VizBin coordinates.')
    coord_df = pd.read_csv(vizbin_coords, sep='\t', names=['contig', 'x', 'y'])
    print('Reading average depth.')
    depth_df = pd.read_csv(depth_file, sep='\t', names=['contig', 'depth'])
    print('Merging data')
    # Merge annotation data and VizBin coords first, keeping only contigs with coords, then merge with depth data
    cont_data_df = annot_df.merge(coord_df, how='right', on='contig').merge(depth_df, how='inner', on='contig').sort_values(by='contig')
    return cont_data_df


def knn_vizbin_coords(contig_info_df, pk):
    # scale/normalize?
    knn_c = KNeighborsClassifier(n_neighbors=pk, weights='distance')
    knn_c.fit(contig_info_df.loc[:, ['x', 'y']], contig_info_df.loc[:, 'contig'])
    # sort distance of neighbouring points
    skn = pd.Series(sorted([row[-1] for row in knn_c.kneighbors()[0]]))
    # calculate running standard deviation between 10 neighbouring points
    sdkn = skn.rolling(10, center=True, min_periods=0).std()
    # find the first jump in distances at the higher end
    try:
        est = sorted([e for i, e in zip(sdkn, skn) if i > sdkn.quantile(.975) and e > skn.median()])[0]
    except IndexError:
        est = None
    return est


def contig_df2cluster_dict(contig_info_df, dbscan_labels):
    # This function is garbage. From pandas df to dict to df. Need to improve asap
    cluster_dict = {}
    contig_info_df = contig_info_df.loc[:, ['contig', 'essential', 'depth', 'x', 'y']]
    contig_info_df['cluster'] = dbscan_labels
    tmp_contig_dict = contig_info_df.fillna('non_essential').set_index('contig').to_dict('index')
    for contig in tmp_contig_dict:
        contig_cluster = shorten_cluster_names(str(tmp_contig_dict.get(contig, {}).get('cluster')))
        contig_essential = tmp_contig_dict.get(contig, {}).get('essential')
        contig_depth = tmp_contig_dict.get(contig, {}).get('depth')
        contig_x = tmp_contig_dict.get(contig, {}).get('x')
        contig_y = tmp_contig_dict.get(contig, {}).get('y')
        if not cluster_dict.get(contig_cluster) and '-1' not in contig_cluster:
            cluster_dict[contig_cluster] = {'essential': [contig_essential],  # .split(','),
                                            'depth': [contig_depth],
                                            'contigs': [contig],
                                            'x': [contig_x],
                                            'y': [contig_y]}
        elif '-1' not in contig_cluster:
            if contig_essential:
                cluster_dict.get(contig_cluster, {}).get('essential').append(contig_essential)  # .extend(contig_essential.split(','))
            cluster_dict.get(contig_cluster, {}).get('depth').append(contig_depth)
            cluster_dict.get(contig_cluster, {}).get('contigs').append(contig)
            cluster_dict.get(contig_cluster, {}).get('x').append(contig_x)
            cluster_dict.get(contig_cluster, {}).get('y').append(contig_y)
    return cluster_dict


def dbscan_cluster(contig_data_df, pk=None, n_jobs=1):
    if not pk:
        pk = int(np.log(contig_data_df['contig'].size))
    # Get reachability distance estimate
    est = knn_vizbin_coords(contig_data_df, pk)
    # Run parallelized dbscan
    df_vizbin_coordinates = contig_data_df.loc[:, ['x', 'y']].to_numpy(dtype=np.float64)
    # print('Running dbscan.')
    with parallel_backend('threading'):
        dbsc = DBSCAN(eps=est, min_samples=pk, n_jobs=n_jobs).fit(df_vizbin_coordinates)
    cluster_labels = dbsc.labels_
    cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)
    return cluster_dict, cluster_labels


def optics_cluster(contig_data_df, min_samples=None, include_depth=False, n_jobs=1):
    if not include_depth:
        dim_df = contig_data_df.loc[:, ['x', 'y']].to_numpy(dtype=np.float64)
    else:
        dim_df = contig_data_df.loc[:, ['x', 'y', 'depth']].to_numpy(dtype=np.float64)
    if not min_samples:
        min_samples = int(np.log(contig_data_df['contig'].size))
    # print('Running OPTICS.')
    clustering = OPTICS(min_samples=min_samples, n_jobs=n_jobs).fit(dim_df)
    cluster_labels = clustering.labels_
    cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)
    return cluster_dict, cluster_labels


def run_initial_scan(contig_data_df, initial_cluster_mode, dbscan_threads, pk=None):
    if initial_cluster_mode == 'DBSCAN' or not initial_cluster_mode:
        # Run parallelized dbscan
        first_clust_dict, labels = dbscan_cluster(contig_data_df, pk=pk, n_jobs=dbscan_threads)

    elif initial_cluster_mode == 'OPTICS':
        # Run OPTICS
        first_clust_dict, labels = optics_cluster(contig_data_df, include_depth=True)
    # Write first scan to file
    first_clust_df = cluster_df_from_dict(first_clust_dict)
    first_clust_df.to_csv('contigs2clusters_initial.tsv', sep='\t', index=False)
    return first_clust_dict, labels


def cluster_df_from_dict(cluster_dict):
    cluster_df = pd.DataFrame()
    cluster_df['cluster'] = [cluster for cluster in cluster_dict if not cluster == '-1']
    for metric in ['contigs', 'essential']:
        metric_list = []
        metric_uniq_list = None
        if metric == 'essential':
            metric_uniq_list = []
        for cluster in cluster_dict:
            essential_list = [i for e in cluster_dict.get(cluster, {}).get(metric) for i in e.split(',')
                              if not e == 'non_essential']
            metric_list.append(len(essential_list))
            if metric == 'essential':
                metric_uniq_list.append(len(set(essential_list)))
        cluster_df[metric] = metric_list
        if metric_uniq_list:
            cluster_df['unique_' + metric] = metric_uniq_list
    cluster_df['completeness'] = round(cluster_df['unique_essential'] / 115, 2)
    cluster_df['purity'] = round(cluster_df['unique_essential'] / cluster_df['essential'], 2)
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
    for cluster in cluster_dict:
        sorted_dict[cluster] = {}
        deps = np.array(cluster_dict[cluster]['depth']).argsort()
        for metric in ['depth', 'contigs', 'essential', 'x', 'y']:
            metric_np = np.array(cluster_dict[cluster][metric])
            sorted_dict[cluster][metric] = metric_np[deps]
    return sorted_dict


def gather_cluster_data(cluster, cluster_dict):
    cluster_essential_genes = [gene for genes in cluster_dict.get(cluster, {}).get('essential')
                               for gene in genes.split(',') if not gene == 'non_essential' ]
    cluster_unique_essential_genes = set(cluster_essential_genes)
    if cluster_essential_genes:
        cluster_purity = round(len(cluster_unique_essential_genes) / len(cluster_essential_genes), 2)
        cluster_completeness = round(len(cluster_unique_essential_genes) / 115, 2)
    else:
        cluster_purity = 0
        cluster_completeness = 0
    cluster_data = [cluster_dict.get(cluster, {}).get('contigs'), cluster_dict.get(cluster, {}).get('essential'),
                    cluster_dict.get(cluster, {}).get('depth'), cluster_dict.get(cluster, {}).get('x'),
                    cluster_dict.get(cluster, {}).get('y'), cluster_purity, cluster_completeness]
    return cluster_data


def create_new_clusters(cluster, intervals, clust_dat):
    new_clusters = {}
    assigned_indices = []
    bin_id = 1
    sum_sub_cluster_indices = 0
    for interval in intervals:
        cluster_id = cluster + '.' + str(bin_id)
        new_clusters[cluster_id] = {'depth': clust_dat[2][interval[0]:interval[1]],
                                    'essential': clust_dat[1][interval[0]:interval[1]],
                                    'contigs': clust_dat[0][interval[0]:interval[1]],
                                    'x': clust_dat[3][interval[0]:interval[1]],
                                    'y': clust_dat[4][interval[0]:interval[1]]}
        bin_id += 1
        # remove indices from leftover list
        assigned_indices.extend(range(interval[0], interval[1]))
        sum_sub_cluster_indices += len(clust_dat[2][interval[0]:interval[1]])
    # Write leftover 'cluster'
    cluster_id = cluster + '.L'
    new_clusters[cluster_id] = {'depth': np.array([clust_dat[2][index] for index, i in enumerate(clust_dat[2].tolist())
                                                   if index not in assigned_indices]),
                                'essential': np.array(
                                    [clust_dat[1][index] for index, i in enumerate(clust_dat[1].tolist())
                                     if index not in assigned_indices]),
                                'contigs': np.array(
                                    [clust_dat[0][index] for index, i in enumerate(clust_dat[0].tolist())
                                     if index not in assigned_indices]),
                                'x': np.array([clust_dat[3][index] for index, i in enumerate(clust_dat[3].tolist())
                                               if index not in assigned_indices]),
                                'y': np.array([clust_dat[4][index] for index, i in enumerate(clust_dat[4].tolist())
                                               if index not in assigned_indices])}
    bin_id += 1
    return new_clusters


def get_sub_clusters(cluster, cluster_dict, threads_for_dbscan, purity_threshold=0.9, completeness_threshold=0.4, alpha=0.001, pk=None, cluster_mode=None, include_depth=False):
    # Create dict with just cluster and sort again, to ensure order by depth is
    cluster_dict = sort_cluster_dict_data_by_depth({cluster: cluster_dict[cluster]})
    # All data needed stored in list with following order:
    # 0:contigs, 1:genes, 2:depth, 3:x, 4:y, 5:purity, 6:completeness
    clust_dat = gather_cluster_data(cluster, cluster_dict)
    if clust_dat[5] < purity_threshold and isinstance(clust_dat[5], float) and clust_dat[6] > completeness_threshold:
        print('Cluster {0} below purity threshold of {1} with {2} and above completeness threshold of {3} with {4}. '
              'Attempting to split.'.format(shorten_cluster_names(cluster), purity_threshold, clust_dat[5], completeness_threshold, clust_dat[6]))
        intervals = [1]
        while len(intervals) == 1 and alpha < 0.05:
            intervals = UniDip(clust_dat[2], alpha=alpha, ntrials=100, mrg_dst=1).run()
            alpha += 0.001
            if len(intervals) > 1:
                print('Found {0} depth sub-clusters in cluster {1} with alpha of {2}'.format(
                    len(intervals), shorten_cluster_names(cluster), round(alpha, 3)))
                new_clusters = create_new_clusters(cluster, intervals, clust_dat)
                return [new_clusters, cluster]
        print('Failed to find depth sub-bclusters for {0}'.format(shorten_cluster_names(cluster)))
        # initialise some stuff
        cluster_contig_df = contig_df_from_cluster_dict(cluster_dict)
        if not pk:
            pk = int(np.log(cluster_contig_df['contig'].size))
            # if pk < 10 < int(np.log(cluster_contig_df['contig'].size)):
            #     pk = 10
        new_clusters_labels = [1]
        if cluster_mode == 'DBSCAN' or not cluster_mode:
            print('Trying with DBSCAN.')
            while cluster_contig_df['contig'].size < pk and pk > 0:
                pk = int(pk * 0.75)
            dbscan_tries = 0
            print('Working on {0}.'.format(shorten_cluster_names(cluster)))
            while len(set(new_clusters_labels)) == 1 and 0 < pk < cluster_contig_df['contig'].size:
                dbscan_tries += 1
                cluster_est = knn_vizbin_coords(cluster_contig_df, pk)
                while not cluster_est and pk > 0:
                    pk = int(pk * 0.75)
                    if pk > 0:
                        cluster_est = knn_vizbin_coords(cluster_contig_df, pk)
                if cluster_est and cluster_contig_df['contig'].size > pk > 0:
                    df_vizbin_coords = cluster_contig_df.loc[:, ['x', 'y']].to_numpy(dtype=np.float64)
                    with parallel_backend('threading'):
                        dbsc = DBSCAN(eps=cluster_est, min_samples=pk, n_jobs=threads_for_dbscan).fit(df_vizbin_coords)
                    new_clusters_labels = dbsc.labels_
                    if len(set(new_clusters_labels)) > 1:
                        new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
                                             enumerate(set(new_clusters_labels))}
                        new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
                        print('Found {0} sub-clusters in cluster {1} with pk of {2}'.format(
                            len(set(new_clusters_labels)), cluster, pk))
                        new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
                        return [new_clusters, cluster]
                    pk = int(pk * 0.75)
                else:
                    print('Failed to find sub-bclusters for {0} with DBSCAN: Could not estimate reachability distance.'.format(shorten_cluster_names(cluster)))
                    return
        elif cluster_mode == "OPTICS":
            print('Trying with OPTICS.')
            # Try with automatic pk first
            min_samples = pk
            optics_tries = 0
            try_interval = 1
            while len(set(new_clusters_labels)) == 1 and 0 < min_samples < cluster_contig_df['contig'].size and optics_tries <= 100:
                if optics_tries == try_interval * 10:
                    try_interval += 1
                    min_samples = int(min_samples / 2)
                    print('Currently at OPTICS try {0} with min samples of {1} for cluster {2}.'.format(optics_tries,
                                                                                                               min_samples,
                                                                                                               shorten_cluster_names(cluster)))
                new_clusters_labels = optics_cluster(cluster_contig_df, min_samples=min_samples, include_depth=include_depth,
                                                     n_jobs=threads_for_dbscan)[1]
                if len(set(new_clusters_labels)) > 1:
                    new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
                                         enumerate(set(new_clusters_labels))}
                    new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
                    print('Found {0} sub-clusters in cluster {1} with min samples of {2}'.format(
                        len(set(new_clusters_labels)), cluster, min_samples))
                    new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
                    return [new_clusters, cluster]
                optics_tries += 1
                min_samples = int(min_samples * 0.75)
            print('Failed to find sub-bclusters for {0} with OPTICS.'.format(shorten_cluster_names(cluster)))
            return
    else:
        if clust_dat[5] == 0:
            print('Could not calculate purity for cluster {0}. Leaving at 0 and skipping.'.format(shorten_cluster_names(cluster)))
            return [{}, cluster]
        elif clust_dat[6] < completeness_threshold:
            print('Cluster {0} below completeness threshold with {1}. Skipping.'.format(shorten_cluster_names(cluster),
                                                                                        clust_dat[6]))
            return [{}, cluster]
        else:
            print('Cluster {0} meets purity threshold of {1} with {2} and has completeness of {3}.'.format(shorten_cluster_names(cluster),
                                                                                                           purity_threshold,
                                                                                                           clust_dat[5],
                                                                                                           clust_dat[6]))


def divide_clusters_by_depth(ds_clstr_dict, threads, min_purity=0.9, min_completeness=0.4, pk=None, cluster_mode=None, include_depth=False):
    min_purity = min_purity / 100
    min_completeness = min_completeness / 100
    dict_cp = ds_clstr_dict.copy()
    n_tries = 0
    clusters_to_process = list(dict_cp.keys())
    # with Parallel(n_jobs=threads) as parallel:
    inner_max_threads = int(threads/len(clusters_to_process))
    if inner_max_threads < 1:
        inner_max_threads = 1
    with parallel_backend("loky", inner_max_num_threads=inner_max_threads):
        while n_tries < 100 and clusters_to_process:
            n_tries += 1
            print('Try: {0}'.format(n_tries))
            sub_clstr_res = Parallel(n_jobs=len(clusters_to_process))(delayed(get_sub_clusters)(str(cluster), dict_cp,
                                                                                                inner_max_threads,
                                                                                                purity_threshold=min_purity,
                                                                                                completeness_threshold=min_completeness,
                                                                                                pk=pk, cluster_mode=cluster_mode,
                                                                                                include_depth=include_depth)
                                                                      for cluster in clusters_to_process)
            for output in sub_clstr_res:
                if output:
                    print('Removing {0}.'.format(output[1]))
                    dict_cp.pop(output[1])
                    print('Adding sub-clusters of {0}.'.format(output[1]))
                    dict_cp.update(output[0])
            clusters_to_process = [i for e in sub_clstr_res if e for i in list(e[0].keys())]
    print('Tries until end: {0}'.format(n_tries))
    return dict_cp


def contig_df_from_cluster_dict(cluster_dict):
    clust_contig_df_rows = []
    clust_contig_df_cols = ['contig', 'essential', 'x', 'y', 'depth', 'cluster', 'purity', 'completeness']
    clust_contig_df_ind = []
    clust_contig_df_ind_init = 0
    for i in cluster_dict:
        cluster_essential_genes = [gene for genes in cluster_dict.get(i, {}).get('essential')
                                   for gene in genes.split(',') if not gene == 'non_essential']
        cluster_unique_essential_genes = set(cluster_essential_genes)
        if cluster_essential_genes:
            cluster_purity = round(len(cluster_unique_essential_genes) / len(cluster_essential_genes), 2)
            cluster_completeness = round(len(cluster_unique_essential_genes) / 115, 2)
        else:
            cluster_purity = 0
            cluster_completeness = 0
        for index, contig in enumerate(cluster_dict[i]['contigs']):
            new_row = [contig, cluster_dict[i]['essential'][index], cluster_dict[i]['x'][index],
                       cluster_dict[i]['y'][index], cluster_dict[i]['depth'][index], i, cluster_purity, cluster_completeness]
            clust_contig_df_rows.append(new_row)
            clust_contig_df_ind.append(clust_contig_df_ind_init)
            clust_contig_df_ind_init += 1
    clust_contig_df = pd.DataFrame(clust_contig_df_rows, clust_contig_df_ind,
                                          clust_contig_df_cols)
    return clust_contig_df


def write_scatterplot(df, file_path, hue):
    palette = sns.color_palette("husl", len(set(hue)))
    sorted_set = []
    for i in hue:
        if not i in sorted_set:
            sorted_set.append(i)
    for index, i in enumerate(sorted_set):
        if i == 'N':
            palette[index] = ('silver')
    scatter_plot = sns.scatterplot(data=df, x="x", y="y", hue=hue, sizes=0.4, s=4, palette=palette)
    scatter_plot.legend(fontsize=3, title="Clusters", title_fontsize=4, ncol=1,
                        bbox_to_anchor=(1.01, 1), borderaxespad=0)
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
        while index+counter+1 <= len(sub_clusters)-1 and sub_clstr == sub_clusters[index+counter+1]:
            counter += 1
        if counter > 0:
            new_name.append(sub_clstr+'x'+str(counter+1))
        else:
            new_name.append(sub_clstr)
    new_name = '.'.join(new_name)
    return new_name


def write_bins(cluster_dict, assembly, min_comp=25, min_pur=90, bin_dir='bins'):
    assembly_dict = load_fasta(assembly)
    # Create bin folder, if it doesnt exist.
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    for cluster in cluster_dict:
        cluster_essential_genes = [gene for genes in cluster_dict.get(cluster, {}).get('essential')
                                   for gene in genes.split(',') if not gene == 'non_essential']
        cluster_unique_ssential_genes = set(cluster_essential_genes)
        if cluster_essential_genes:
            cluster_purity = round(len(cluster_unique_ssential_genes) / len(cluster_essential_genes) * 100, 1)
            cluster_completeness = round(len(cluster_unique_ssential_genes) / 115 * 100, 1)
        else:
            cluster_purity = 0
            cluster_completeness = 0
        if cluster_purity >= min_pur and cluster_completeness >= min_comp:
            new_cluster_name = shorten_cluster_names(cluster)
            bin_name = '_'.join(['binny'] + ['C'+str(cluster_completeness)] + ['P'+str(cluster_purity)] + [new_cluster_name])
            bin_file_name = bin_name + '.fasta'
            bin_out_path = bin_dir / bin_file_name
            with open(bin_out_path, 'w') as out_file:
                for contig in cluster_dict[cluster]['contigs']:
                    out_file.write('>' + contig + '\n' + assembly_dict.get(contig) + '\n')
