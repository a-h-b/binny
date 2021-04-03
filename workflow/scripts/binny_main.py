#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:50:35 2020

@author: oskar.hickl
"""

cluster_dir = snakemake.input['outdir']
mg_depth_file = snakemake.input['mgdepth']
assembly = snakemake.input['assembly']
coords_file = snakemake.input['vizbin']
annot_file = snakemake.input['gff']
functions = snakemake.params['py_functions']
<<<<<<< Updated upstream
pk = snakemake.config['binning']['binny']['pk']
all_outputs = snakemake.output
min_purity = snakemake.params['purity']
min_completeness = snakemake.params['completeness']
=======
min_purity = int(snakemake.params['purity'])
min_completeness = int(snakemake.params['completeness'])
kmers = snakemake.params['kmers']
min_contig_length = int(snakemake.params['cutoff'])
>>>>>>> Stashed changes
threads = snakemake.threads
log = snakemake.log[0]

import sys
sys.path.append(functions)
from binny_functions import *

sys.stdout = open(log, 'w')

<<<<<<< Updated upstream
=======
# Load assembly and mask rRNAs and CRISPR arrays
contig_list = [[key] + [val] for key, val in load_fasta(assembly).items() if len(val) >= min_contig_length]
contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)

# Get length normalized k-mer frequencies.
kmer_sizes = [int(kmer) for kmer in kmers.split(',')]
print(kmer_sizes)
start = timer()
kfreq_array = get_contig_kmer_matrix2(contig_list, kmer_sizes, threads)
end = timer()
print('K-mer frequency matrix created in {0}s.'.format(end - start))

# Standardize k-mer freq data
X = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:]])
X_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:]]
scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)

# Manifold learning and dimension reduction.
print('Running manifold learning and dimension reduction')
n_dim = 3
perp = get_perp(len(X_contigs))
# print('Perplexity: {0}.'.format(perp))
X_embedded = umap.UMAP(n_neighbors=perp, min_dist=0.0, n_components=n_dim, random_state=0, densmap=False, verbose=True).fit_transform(X_scaled)

# Create 2/3D coordinate df.
if n_dim == 3:
    coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['x', 'y', 'z'])
    coord_df['contig'] = X_contigs
    coord_df = coord_df[['contig', 'x', 'y', 'z']]
    coord_df.to_csv('intermediary/contig_coordinates.tsv', sep='\t', index=False, header=False)
    coords_file = 'intermediary/contig_coordinates.tsv'
elif n_dim == 2:
    coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['x', 'y'])
    coord_df['contig'] = X_contigs
    coord_df = coord_df[['contig', 'x', 'y']]
    coord_df.to_csv('intermediary/contig_coordinates.tsv', sep='\t', index=False, header=False)
    coords_file = 'intermediary/contig_coordinates.tsv'

>>>>>>> Stashed changes
# Load data
contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, coords_file)

# Write contig data to file
contig_data_df.to_csv('contig_data.tsv', sep='\t', index=False)

# Run initial clustering
first_clust_dict, labels = run_initial_scan(contig_data_df, 'DBSCAN', threads)

# Plot initial clustering
write_scatterplot(contig_data_df, 'initial_scatter_plot.pdf', labels)

# Find sub-clusters
print('Attempting to divide clusters by depth and/or sub-clusters by subsequent dbscan/optics runs.')
final_clust_dict = divide_clusters_by_depth(first_clust_dict, threads, int(min_purity), int(min_completeness),
                                            cluster_mode='OPTICS',  include_depth=True)
final_clust_df = cluster_df_from_dict(final_clust_dict)
final_clust_df.to_csv('contigs2clusters_final.tsv', sep='\t', index=False)
final_clust_contig_df = contig_df_from_cluster_dict(final_clust_dict)

# Plot final clustering
conditions = [(final_clust_contig_df['purity'] >= min_purity) & (final_clust_contig_df['completeness'] >= min_completeness),
              (final_clust_contig_df['purity'] < min_purity) | (final_clust_contig_df['completeness'] < min_completeness)]
values = [final_clust_contig_df['cluster'], 'N']
final_clust_contig_df['above_thresh'] = np.select(conditions, values)

final_clust_contig_df = contig_data_df.merge(final_clust_contig_df, how='outer', on='contig', suffixes=(None, '_y'))
final_clust_contig_df['above_thresh'] = final_clust_contig_df['above_thresh'].fillna('N')

write_scatterplot(final_clust_contig_df, 'final_scatter_plot.pdf', final_clust_contig_df['above_thresh'])

write_bins(final_clust_dict, assembly, min_comp=int(min_completeness), min_pur=int(min_purity), bin_dir='bins')

print('Run finished.')
sys.stdout.close()
