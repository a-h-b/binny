#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:50:35 2021

@author: oskar.hickl
"""

import logging
import sys
import os

binny_out = snakemake.params['binny_out']
sample = snakemake.params['sample']
mg_depth_file = snakemake.input['mgdepth']
assembly = snakemake.input['assembly']
annot_file = snakemake.params['gff']
raw_annot = snakemake.input['raw_gff']
tigrfam2pfam_file = snakemake.params['t2p']
taxon_marker_set_file = snakemake.params['marker_sets']
prokka_checkm_marker_hmm_out = snakemake.input['hmm_markers']
functions = snakemake.params['py_functions']
min_purity = float(snakemake.params['purity'])
min_completeness = float(snakemake.params['min_completeness'])
starting_completeness = float(snakemake.params['start_completeness'])
kmers = snakemake.params['kmers']
nx_val = int(snakemake.params['nx_val'])
min_contig_length = int(snakemake.params['min_cutoff'])
max_contig_length = int(snakemake.params['max_cutoff'])
min_contig_length_marker = int(snakemake.params['min_cutoff_marker'])
max_contig_length_marker = int(snakemake.params['max_cutoff_marker'])
max_contig_threshold = float(snakemake.params['max_n_contigs'])
max_embedding_tries = int(snakemake.params['max_embedding_tries'])
include_depth_initial = eval(snakemake.params['include_depth_initial'])
include_depth_main = eval(snakemake.params['include_depth_main'])
hdbscan_epsilon_range = [float(epsilon) for epsilon in snakemake.params['hdbscan_epsilon_range'].split(',')]
hdbscan_min_samples_range = [int(min_sample) for min_sample in snakemake.params['hdbscan_min_samples_range'].split(',')]
dist_metric = snakemake.params['distance_metric']

intermediary_file_dir = 'intermediary'

threads = snakemake.threads
log = snakemake.log[0]

n_dim = 2

sys.path.append(functions)
from binny_functions import *

# To achieve reproducible results with HDBSCAN and ensure same seed, because other tools that accept seed arguments,
# might mess with the global numpy seed
np.random.seed(0)

logging.basicConfig(filename=log, level=logging.INFO, format='%(asctime)s - %(message)s',  # logging.INFO
                    datefmt='%d/%m/%Y %I:%M:%S %p', filemode='w')

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

logging.info('Starting Binny run for sample {0}.'.format(sample))

# Load TIGRFAMs to PFAMs conversion table.
tigrfam2pfam_data = tigrfam2pfam_dict(tigrfam2pfam_file)

# Merge Prokka gff with marker set data, load annotation df, and load assembly.
checkm_hmmer_search2prokka_gff(prokka_checkm_marker_hmm_out, raw_annot, gff_out_path=os.path.join(binny_out,
                                                                                                  'intermediary'))
annot_df, annot_dict = gff2ess_gene_df(annot_file, target_attribute='checkm_marker', get_dict=True)
assembly_dict = load_fasta(assembly)

# Build marker set graph db.
taxon_marker_sets = load_checkm_markers(taxon_marker_set_file)

# Look for complete genomes on single contigs
all_good_bins = {}
logging.info('Looking for single contig bins.')
single_contig_bins = get_single_contig_bins(annot_df, all_good_bins, n_dim, taxon_marker_sets, tigrfam2pfam_data,
                                            threads)
logging.info('Found {0} single contig bins.'.format(len(single_contig_bins)))

logging.info('Getting assembly dict without scMAG.'.format(len(single_contig_bins)))
# assembly_dict_wo_scmags = {key: val for key, val in assembly_dict.items() if key in
#                            set(assembly_dict.keys()).difference(set(single_contig_bins))}

logging.info(f'Calculating N{nx_val}'.format(len(single_contig_bins)))
nx = calc_assembly_nx(assembly_dict, single_contig_bins, nx_val)
nx2 = calc_assembly_nx(assembly_dict, [], nx_val)
logging.info(f'N{nx_val} is {nx}, with scMAGs would be {nx2}.'.format(len(single_contig_bins)))
min_contig_length = min(max(nx, min_contig_length), max_contig_length)
min_contig_length_marker = min(max(int(nx / 3), min_contig_length_marker), max_contig_length_marker)

# Load assembly and mask rRNAs and CRISPR arrays
contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if (len(seq) >= min_contig_length
                                                                             or (annot_dict.get(contig)
                                                                             and len(seq) >= min_contig_length_marker))
                                                                            and contig not in single_contig_bins]

logging.info('{0} contigs match length threshold of {1}bp or contain marker genes and'
             ' have a size of at least {2}bp'.format(len(contig_list), min_contig_length, min_contig_length_marker))

contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)  # Disabled for v016

# Get length normalized k-mer frequencies.
kmer_sizes = [int(kmer) for kmer in kmers.split(',')]

start = timer()
kfreq_array = get_contig_kmer_matrix(contig_list, kmer_sizes, threads)
end = timer()
logging.info('K-mer frequency matrix created in {0}s.'.format(int(end - start)))

# Make array, removing fully masked sequences with no counts and standardize k-mer freq data
x = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0])
x_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0]

main_contig_data_dict = {cont: seq for cont, seq in zip(x_contigs, x)}

# Load depth data
depth_dict = load_depth_dict(mg_depth_file)

# Run iterative dimension reduction, manifold learning, cluster detection and assessment.
all_good_bins, contig_data_df_org, min_purity = iterative_embedding(x_contigs, depth_dict, all_good_bins, starting_completeness,
                                                        min_purity, min_completeness, threads, n_dim, annot_file,
                                                        mg_depth_file, single_contig_bins, taxon_marker_sets,
                                                        tigrfam2pfam_data, main_contig_data_dict, assembly_dict,
                                                        max_contig_threshold, min_contig_length_marker,
                                                        include_depth_initial, max_embedding_tries,
                                                        include_depth_main, hdbscan_epsilon_range,
                                                        hdbscan_min_samples_range, dist_metric,
                                                        contigs2clusters_out_path=os.path.join(binny_out,
                                                                                               'intermediary'))

all_contigs = []
for bin in all_good_bins:
    all_contigs.extend(all_good_bins[bin]['contigs'])
if len(all_contigs) != len(set(all_contigs)):
    logging.warning('WARNING: {0} duplicate contigs in bins found!'.format(len(all_contigs) - len(set(all_contigs))))

# Write bin fastas.
write_bins(all_good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity),
           bin_dir=os.path.join(binny_out, 'bins'))

bin_dict = {contig: bin for bin, bin_data in all_good_bins.items() for contig in bin_data['contigs']}

all_cont_data_dict = {}

for contig, k_freq in main_contig_data_dict.items():
    all_cont_data_dict[contig] = {'bin': bin_dict.get(contig, 'N'),
                                  'k-mer_freqs': ';'.join([str(k) for k in list(k_freq)]),
                                  'depths': ';'.join([str(d) for d in list(depth_dict.get(contig))])}

all_cont_data_dict.update({contig: {'bin': contig,
                                  'k-mer_freqs': '',
                                  'depths': ';'.join([str(d) for d in list(depth_dict.get(contig))])}
                           for contig in single_contig_bins})

compression_opts = dict(method='gzip')

contig_data_df = pd.DataFrame.from_dict(all_cont_data_dict, orient='index', columns=['bin', 'k-mer_freqs', 'depths'])

contig_data_df.to_csv(os.path.join(binny_out, 'contig_data.tsv.gz'), header=True, index=True, index_label='contig',
                      chunksize=1000, compression=compression_opts, sep='\t')

logging.info('Run finished.')
