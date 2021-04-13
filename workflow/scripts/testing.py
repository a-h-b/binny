binny_functions_path = '/Users/oskar.hickl/local_tools/binny_fork/binny/workflow/scripts/'

import os
import sys

sys.path.append(binny_functions_path)
from binny_functions import *

########################################################################################################################
os.chdir('/Users/oskar.hickl/binny_bench/binner_data/Binny/binny_outputs/2017.12.04_18.45.54_sample_5')
########################################################################################################################
binnydir = 'intermediary/'
mg_depth_file = 'intermediary/assembly.contig_depth.txt'
assembly = 'assembly.fa'
# annot_file = 'intermediary/annotation_CDS_RNA_hmms.gff'
# annot_file = 'intermediary/annotation_CDS_RNA_hmms_checkm_test2.gff'  # for smaple 5
annot_file = 'intermediary/annotation_CDS_RNA_hmms_checkm.gff'
raw_annot = 'intermediary/annotation.filt.gff'
# prokka_checkm_marker_hmm_out = 'intermediary/prokka.faa.checkm.v3.tblout.hmmscan' # for smaple 5
# prokka_checkm_marker_hmm_out = 'intermediary/prokka.faa.checkm.hmmscan'
prokka_checkm_marker_hmm_out = 'intermediary/prokka.faa.checkm.domtblout.v4.hmmscan'
# annot_file_checkm = 'intermediary/annotation_CDS_RNA_hmms_checkm_test.gff'
tigrfam2pfam_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/pfam/tigrfam2pfam.tsv'
taxon_marker_set_file = '/Users/oskar.hickl/Downloads/checkm_data_2015_01_16/taxon_marker_sets_lineage_sorted.tsv'

starting_completeness = 90  # 90
min_completeness = 70  # 70
min_purity = 90  # 90
threads = 10
n_dim = 2

tigrfam2pfam_data = tigrfam2pfam_dict(tigrfam2pfam_file)

# annot_df, annot_dict = gff2ess_gene_df(annot_file, get_dict=True)
checkm_hmmer_search2prokka_gff_v2(prokka_checkm_marker_hmm_out, raw_annot, tigrfam2pfam_data)
annot_df, annot_dict = gff2ess_gene_df(annot_file, target_attribute='checkm_marker', get_dict=True)
assembly_dict = load_fasta(assembly)

taxon_marker_sets = load_checkm_markers(taxon_marker_set_file)

# Look for complete genomes on single contigs
all_good_bins = {}
print('Looking for single contig bins.')
single_contig_bins = get_single_contig_bins(annot_df, all_good_bins, n_dim, taxon_marker_sets, tigrfam2pfam_data, threads=10)
print('Found {0} single contig bins.'.format(len(single_contig_bins)))
print(len(all_good_bins.keys()))

min_cont_size = 1000

# Load assembly and mask rRNAs and CRISPR arrays
# contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if (len(seq) >= min_cont_size
#                                                                        or (annot_dict.get(contig) and len(seq) >= min_cont_size / 2))
#                                                                         and contig not in single_contig_bins]
contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if (len(seq) >= min_cont_size or annot_dict.get(contig))
                                                                        and contig not in single_contig_bins]
print(len(contig_list))
contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)

# Get length normalized k-mer frequencies.
print('Counting k-mers.')
kmer_sizes = [2, 3, 4]  # list(range(6, 7))
start = timer()
kfreq_array = get_contig_kmer_matrix2(contig_list, kmer_sizes, threads)
end = timer()
print('K-mer frequency matrix created in {0}s.'.format(end - start))

# Make array, removing fully masked sequences with no counts and standardize k-mer freq data
X = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0])
X_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0]

main_contig_data_dict = {cont: seq for cont, seq in zip(X_contigs, X)}

########################################################################################################################
gs = 'gs_cont2gen.tsv'
cont2gen_gs_dict = {}
with open(gs, 'r') as f:
    for line in f:
        line = line.strip('\n \t').split('\t')
        cont2gen_gs_dict[line[0]] = line[1]
########################################################################################################################

depth_dict = load_depth_dict(mg_depth_file)
########################################################################################################################
embedding_tries = 1
early_exag = 100
internal_completeness = starting_completeness
########################################################################################################################
while embedding_tries <= 100:
    round_X_contigs = []
    round_X = []
    print('Running embedding iteration {0}.'.format(embedding_tries))
    if embedding_tries == 1:
        round_X_contigs = X_contigs
        round_X = X
    else:
        round_X = [main_contig_data_dict.get(cont) for cont in round_leftovers['contig'].tolist()]
        round_X_contigs = [cont for cont in round_leftovers['contig'].tolist()]
        # for index, i in enumerate(X_contigs):
        #     if i in round_leftovers['contig'].tolist():
        #         round_X_contigs.append(X_contigs[index])
        #         round_X.append(X[index])
    if len(round_X_contigs) != len(round_X):
        print('fuuuuuuuuck')
        raise Exception
    # print(len(round_X_contigs))
    # round_X_contigs = round_X
    # Replace zeroes for clr
    # round_X = np.array([np.append(np.array(depth_dict[contig]), kmers) for contig, kmers in zip(round_X_contigs, round_X)])
    round_X = multiplicative_replacement(round_X)
    # Clr transform
    X_scaled = clr(round_X)
    X_scaled = np.array([np.append(np.array(depth_dict[contig]), kmers) for contig, kmers in zip(round_X_contigs, X_scaled)])

    # os.chdir('/Users/oskar.hickl/local_tools/EMBEDR')
    # X_scaled = np.loadtxt('/Users/oskar.hickl/local_tools/EMBEDR/Data/mnist2500_X.txt')

    # Manifold learning and dimension reduction.
    print('Running manifold learning and dimension reduction')
    # perp = int(get_perp(len(round_X_contigs)))  # / embedding_tries) # - embedding_tries
    # perp = 3000

    n_comp = None
    sum_var_exp = 0
    n_tries = 0
    if len(X_scaled[0]-1) < 25:
        n_comp = len(X_scaled[0]-1)
    else:
        n_comp = 25

    pca = PCA(n_components=n_comp, random_state=0)
    transformer = pca.fit(X_scaled)
    sum_var_exp = sum(pca.explained_variance_ratio_)
    # print(n_comp, sum(pca.explained_variance_ratio_), pca.explained_variance_ratio_)

    while sum_var_exp <= 0.75 and n_tries <= 100 and len(pca.explained_variance_ratio_) <= 75:
        pca = PCA(n_components=n_comp)
        transformer = pca.fit(X_scaled)
        sum_var_exp = sum(pca.explained_variance_ratio_)
        # print(pca.explained_variance_ratio_)
        # print(pca.singular_values_)
        n_comp += 5
        n_tries += 1
    print(n_comp, round(sum(pca.explained_variance_ratio_), 3))  #, pca.explained_variance_ratio_)
    X_pca = transformer.transform(X_scaled)

    if len(round_X_contigs) >= 1e5:
        preplexities = [4, 10, 100, 1000, 5000]
    else:
        # preplexities = [4, 10, 100, 1000]
        preplexities = [4, 10, 100, 500]  # for sample 0, so it fits in ram

    affinities_multiscale_mixture = affinity.Multiscale(X_pca, perplexities=preplexities, metric="manhattan", n_jobs=threads, random_state=0, verbose=True)  # [perp, perp*10, perp*100] #, 10000 # [2, 10, 100, 1000]
    init = initialization.pca(X_pca, random_state=0, verbose=True)
    embedding = TSNEEmbedding(init, affinities_multiscale_mixture, negative_gradient_method="fft", n_jobs=threads, random_state=0, verbose=True)
    embedding1 = embedding.optimize(n_iter=250, exaggeration=early_exag, momentum=0.5, n_jobs=threads, verbose=True)
    embedding2 = embedding1.optimize(n_iter=750, exaggeration=1, momentum=0.8, n_jobs=threads, verbose=True)
    embedding_multiscale = embedding2.view(np.ndarray)
    data_Y = [embedding_multiscale]


    # Plot gs
    # print('Plotting with gold standard genomes')
    # dim_range = [i + 1 for i in range(n_dim)]
    # gs_coord_df = pd.DataFrame(data=data_Y[0], index=None, columns=['dim' + str(i) for i in dim_range])
    # gs_coord_df['contig'] = round_X_contigs
    # gs_coord_df['cluster'] = [cont2gen_gs_dict [contig] for contig in round_X_contigs]
    # gs_coord_df = gs_coord_df[['contig'] + ['dim' + str(i) for i in dim_range] + ['cluster']]
    # write_scatterplot(gs_coord_df, gs_coord_df['cluster'])

    round_X_bad_contigs = None
    round_X_good_contigs = None
    # print(embedr_obj.tSNE_params())

    # #TEST
    # round_X_contigs = X_contigs
    X_embedded = data_Y[0]

    # Create coordinate df.
    dim_range = [i+1 for i in range(n_dim)]
    coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['dim'+str(i) for i in dim_range])
    coord_df['contig'] = round_X_contigs  #_good
    # Reorder
    coord_df = coord_df[['contig']+['dim'+str(i) for i in dim_range]]
    coords_file = 'intermediary/contig_coordinates_dynamic_marker_sets_v06.tsv'
    coord_df.to_csv(coords_file, sep='\t', index=False, header=False)

    # Load data
    contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, coords_file, dims=n_dim)

    for i in single_contig_bins:
        if i in coord_df['contig'].tolist():
            print('TEST')
            raise Exception


    # Write contig data to file
    if embedding_tries == 1:
        contig_data_df_org = contig_data_df.copy()
        contig_data_df_org.to_csv('intermediary/contig_data_original_dynamic_marker_sets_v06.tsv', sep='\t', index=False)
    contig_data_df.to_csv('intermediary/contig_data_dynamic_marker_sets_v06.tsv', sep='\t', index=False)

    # Find bins
    good_bins, final_init_clust_dict = binny_iterate(contig_data_df, threads, taxon_marker_sets, tigrfam2pfam_data,
                                                     min_purity, internal_completeness, 1, embedding_iteration=embedding_tries)

    if not list(good_bins.keys()) and internal_completeness > min_completeness:  # 60
        # if min_completeness < 90:
        #     min_purity = 95
        internal_completeness = internal_completeness - 5
        early_exag = early_exag * 1.1  # * 1.1
        print('Could not find good bins. Lowering completeness threshold to {0} and increasing t-SNE early exaggeration to {1}.'.format(internal_completeness, early_exag))
    elif not list(good_bins.keys()) and not list(round_good_bins.keys()):
        if embedding_tries > 1:
            print('Found no good bins two times in a row')
            break
    elif not list(good_bins.keys()) and internal_completeness < min_completeness:  # 60
        if embedding_tries > 1:
            print('Reached min completeness and found no more bins. Exiting embedding iteration')
            break
    else:
        early_exag = early_exag * 10

    round_good_bins = good_bins

    # if not list(good_bins.keys()):
    #     min_purity = 98
    #     min_completeness = 10
    #     good_bins, final_init_clust_dict = binny_iterate(contig_data_df, threads, min_purity=min_purity, min_completeness=min_completeness,
    #                                                      max_iterations=1, embedding_iteration=embedding_tries)

    print('Good bins this iteration:{0}.'.format(len(good_bins.keys())))
    # Write round table
    round_clust_dict = {**good_bins, **final_init_clust_dict}
    round_clust_df = cluster_df_from_dict(round_clust_dict)
    round_clust_df.to_csv('round_contigs2clusters_dynamic_marker_sets_v06.tsv', sep='\t', index=False)

    round_clust_contig_df = contig_df_from_cluster_dict(round_clust_dict)

    conditions = [(round_clust_contig_df['purity'] >= min_purity / 100) & (
                round_clust_contig_df['completeness'] >= internal_completeness / 100),
                  (round_clust_contig_df['purity'] < min_purity / 100) | (
                              round_clust_contig_df['completeness'] < internal_completeness / 100)]
    values = [round_clust_contig_df['cluster'], 'N']
    round_clust_contig_df['above_thresh'] = np.select(conditions, values)
    round_clust_contig_df = contig_data_df.merge(round_clust_contig_df, how='outer', on='contig', suffixes=(None, '_y'))
    round_clust_contig_df = round_clust_contig_df[round_clust_contig_df.columns.drop(list(round_clust_contig_df.filter(regex='_y')))]
    # print(round_clust_contig_df.columns)
    round_clust_contig_df['above_thresh'] = round_clust_contig_df['above_thresh'].fillna('N')

    round_leftovers = round_clust_contig_df[round_clust_contig_df['above_thresh'] == 'N']

    if round_X_bad_contigs:
        round_leftovers = bad_contig_data_df.merge(round_leftovers, how='outer', on='contig', suffixes=(None, '_y'))
        round_leftovers = round_leftovers[
            round_leftovers.columns.drop(list(round_leftovers.filter(regex='_y')))]

    if len(set(round_leftovers['contig'].tolist())) != len(round_leftovers['contig'].tolist()):
        print('test2')
        print(len(set(round_leftovers['contig'].tolist())), len(round_leftovers['contig'].tolist()))
        raise Exception

    if list(good_bins.keys()):
        good_bin_contig_df = contig_df_from_cluster_dict(good_bins)
        round_leftovers = round_leftovers[~round_leftovers['contig'].isin(good_bin_contig_df['contig'].tolist())]

    # if embedding_tries > 1:
    #     if round_leftovers.equals(round_clust_contig_df[round_clust_contig_df['above_thresh'] == 'N']):
    #         print('test3')
    #         break

    all_good_bins.update(good_bins)
    embedding_tries += 1

    all_contigs = []
    for bin in all_good_bins:
        all_contigs.extend(all_good_bins[bin]['contigs'])
    if len(all_contigs) != len(set(all_contigs)):
        print('WARNING: {0} duplicate contigs in bins found!'.format(len(all_contigs) - len(set(all_contigs))))
        break
    print('Good bins so far:{0}.'.format(len(all_good_bins.keys())))
########################################################################################################################

final_init_clust_dict, labels = run_initial_scan(contig_data_df_org, 'HDBSCAN', threads, include_depth=False)
for cluster in final_init_clust_dict:
    if not final_init_clust_dict[cluster].get('purity'):
        final_init_clust_dict[cluster]['purity'] = 0
        final_init_clust_dict[cluster]['completeness'] = 0
        final_init_clust_dict[cluster]['taxon'] = 'none'

# Write final table
final_clust_dict = {**all_good_bins, **final_init_clust_dict}
# final_clust_dict = {**good_bins, **final_init_clust_dict}
final_clust_df = cluster_df_from_dict(final_clust_dict)
final_clust_df.to_csv('final_contigs2clusters_dynamic_marker_sets_v06.tsv', sep='\t', index=False)

final_clust_contig_df = contig_df_from_cluster_dict(final_clust_dict)

# Plot final clustering
conditions = [(final_clust_contig_df['purity'] >= min_purity / 100) & (
            final_clust_contig_df['completeness'] >= min_completeness / 100),
              (final_clust_contig_df['purity'] < min_purity / 100) | (
                          final_clust_contig_df['completeness'] < min_completeness / 100)]
values = [final_clust_contig_df['cluster'], 'N']
final_clust_contig_df['above_thresh'] = np.select(conditions, values)
final_clust_contig_df = contig_data_df.merge(final_clust_contig_df, how='outer', on='contig', suffixes=(None, '_y'))
final_clust_contig_df = final_clust_contig_df[final_clust_contig_df.columns.drop(list(final_clust_contig_df.filter(regex='_y')))]
# print(final_clust_contig_df.columns)
final_clust_contig_df['above_thresh'] = final_clust_contig_df['above_thresh'].fillna('N')

final_clust_contig_df['gs_genome'] = [cont2gen_gs_dict[contig] for contig in final_clust_contig_df['contig'].tolist()]

########################################################################################################################
gs_stats = gold_standard_stats(cont2gen_gs_dict, assembly_dict)


def asses_bins(final_df, gs_stats):
    final_df_no_n = final_df[final_df['above_thresh'] != 'N']
    final_bins_dict = {}
    for cluster in set(final_df_no_n['cluster'].tolist()):
        cluster_df = final_df[final_df['cluster'] == cluster]
        final_bins_dict[cluster] = {}
        for genome in set(cluster_df['gs_genome'].tolist()):
            genome_size_bp = sum([len(assembly_dict[contig]) for contig in cluster_df['contig'][cluster_df['gs_genome']
                                                                                                == genome].tolist()])
            genome_size_seq = len(cluster_df['contig'][cluster_df['gs_genome'] == genome].tolist())
            final_bins_dict[cluster][genome] = [genome_size_bp, genome_size_seq]
        genome_sizes_bp = [[genome, data[0]] for genome, data in final_bins_dict[cluster].items()]
        genome_sizes_bp.sort(key=lambda x: x[1], reverse=True)
        largest_genome = genome_sizes_bp[0][0]

        final_bins_dict[cluster]['completeness_bp'] = round(final_bins_dict[cluster][largest_genome][0] / gs_stats[largest_genome][0], 3)
        final_bins_dict[cluster]['completeness_seq'] = round(final_bins_dict[cluster][largest_genome][1] / gs_stats[largest_genome][1], 3)

        final_bins_dict[cluster]['purity_bp'] = round(final_bins_dict[cluster][largest_genome][0] / sum([data[0] for genome, data in final_bins_dict[cluster].items() if isinstance(data, list)]), 3)
        final_bins_dict[cluster]['purity_seq'] = round(final_bins_dict[cluster][largest_genome][1] / sum([data[1] for genome, data in final_bins_dict[cluster].items() if isinstance(data, list)]), 3)
        print(cluster, 'Purity[bp]:', final_bins_dict[cluster]['purity_bp'], 'Completeness[bp]:', final_bins_dict[cluster]['completeness_bp'],
              'Purity[seq]:', final_bins_dict[cluster]['purity_seq'], 'Completeness[seq]:', final_bins_dict[cluster]['completeness_seq'])
    return final_bins_dict


bin_stats = asses_bins(final_clust_contig_df, gs_stats)
########################################################################################################################

# if n_dim <= 3:
#     write_scatterplot(final_clust_contig_df, 'final_scatter_plot_dynamic_marker_sets_v06.pdf',
#                       final_clust_contig_df['above_thresh'])

# Check if something went wrong and contigs were duplicated.
all_contigs = []
for bin in all_good_bins:
    all_contigs.extend(all_good_bins[bin]['contigs'])
if len(all_contigs) != len(set(all_contigs)):
    print('WARNING: {0} duplicate contigs in bins found!'.format(len(all_contigs) - len(set(all_contigs))))

# Write bin fastas.
write_bins(all_good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity), bin_dir='bins_dynamic_marker_sets_v06')

print('Run finished.')
