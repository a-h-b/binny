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
pk = snakemake.config['binning']['binny']['pk']
all_outputs = snakemake.output
min_purity = snakemake.params['purity']
min_completeness = snakemake.params['completeness']
threads = snakemake.threads
log = snakemake.log[0]

import sys
sys.path.append(functions)
from binny_functions import *

sys.stdout = open(log, 'w')

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
