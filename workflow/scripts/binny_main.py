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
functions = snakemake.params['plot_functions']
pk = snakemake.config['binning']['binny']['pk']
all_outputs = snakemake.output
min_purity = snakemake.params['purity']
min_completeness = snakemake.params['completeness']
threads = snakemake.threads

import sys
sys.path.append(functions)
from binny_functions import *

# Load data
contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, coords_file)

# Write contig data to file
contig_data_df.to_csv('contig_data.tsv', sep='\t', index=False)

# Run parallelized dbscan
first_clust_dict, labels = dbscan_cluster(contig_data_df, pk, threads)

# Plot first clustering
write_scatterplot(contig_data_df, 'dbscan_scatter_plot_pk{0}.pdf'.format(pk), labels)

# Write first scan to file
first_clust_df = cluster_df_from_dict(first_clust_dict)
first_clust_df.to_csv('contigs2clusters_initial_pk{0}.tsv'.format(pk), sep='\t', index=False)

# first_clust_dict_sorted = sort_cluster_dict_data_by_depth(first_clust_dict)

# Find sub-clusters
print('Attempting to divide dbscan clusters by depth and/or sub-clusters by subsequent dbscan runs.')
final_clust_dict = divide_clusters_by_depth(first_clust_dict, threads)
final_clust_df = cluster_df_from_dict(final_clust_dict)
final_clust_df.to_csv('contigs2clusters_final_pk{0}.tsv'.format(pk), sep='\t', index=False)

final_clust_contig_df = contig_df_from_cluster_dict(final_clust_dict)

# Plot second clustering
write_scatterplot(final_clust_contig_df, 'final_scatter_plot_pk{0}.pdf'.format(pk), final_clust_contig_df['cluster'])

write_bins(final_clust_dict, assembly, min_comp=int(min_completeness), min_pur=int(min_purity))
write_bins(final_clust_dict, assembly, 0, 0, 'bins_no_filter')

print('Run finished.')
