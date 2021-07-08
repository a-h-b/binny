#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:50:35 2021

@author: oskar.hickl
"""

sample = snakemake.params['sample']
# cluster_dir = snakemake.input['outdir']
mg_depth_file = snakemake.input['mgdepth']
assembly = snakemake.input['assembly']
# coords_file = snakemake.input['vizbin']
annot_file = snakemake.params['gff']
raw_annot = snakemake.input['raw_gff']
# hmm_out = snakemake.params['hmm_markers']
tigrfam2pfam_file = snakemake.params['t2p']
taxon_marker_set_file = snakemake.params['marker_sets']
prokka_checkm_marker_hmm_out = snakemake.input['hmm_markers']
functions = snakemake.params['py_functions']
# pk = snakemake.config['binning']['binny']['pk']
# all_outputs = snakemake.output
min_purity = int(snakemake.params['purity'])
min_completeness = int(snakemake.params['completeness'])
kmers = snakemake.params['kmers']
min_contig_length = int(snakemake.params['cutoff'])

threads = snakemake.threads
log = snakemake.log[0]

starting_completeness = 90
min_marker_cont_length = 0
n_dim = 2

import sys
import logging
import pickle
sys.path.append(functions)
from binny_functions import *

# sys.stdout = open(log, 'w')

logging.basicConfig(filename=log, level=logging.DEBUG, format='%(asctime)s - %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p',
                    filemode='w')

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

logging.info('Starting Binny run for sample {0}.'.format(sample))

# Load TIGRFAMs to PFAMs conversion table.
tigrfam2pfam_data = tigrfam2pfam_dict(tigrfam2pfam_file)

# Merge Prokka gff with marker set data, load annotation df, and load assembly.
checkm_hmmer_search2prokka_gff_v2(prokka_checkm_marker_hmm_out, raw_annot, tigrfam2pfam_data)
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

# Load assembly and mask rRNAs and CRISPR arrays
contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if (len(seq) >= min_contig_length
                                                                       or (annot_dict.get(contig) and len(seq) >= min_marker_cont_length))
                                                                        and contig not in single_contig_bins]
logging.info('{0} contigs match length threshold of {1} or contain marker genes and will be used for binning'.format(len(contig_list),
                                                                                                              min_contig_length))
contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)

# Get length normalized k-mer frequencies.
kmer_sizes = [int(kmer) for kmer in kmers.split(',')]
# print('Using k-mer sizes {0}'.format(', '.join(kmer_sizes)))
start = timer()
kfreq_array = get_contig_kmer_matrix2(contig_list, kmer_sizes, threads)
end = timer()
logging.info('K-mer frequency matrix created in {0}s.'.format(int(end - start)))

# Make array, removing fully masked sequences with no counts and standardize k-mer freq data
x = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0])
x_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0]

main_contig_data_dict = {cont: seq for cont, seq in zip(x_contigs, x)}

# Load depth data
depth_dict = load_depth_dict(mg_depth_file)

# Run iterative dimension reduction, manifold learning, cluster detection and assessment.
all_good_bins, contig_data_df_org = iterative_embedding(x, x_contigs, depth_dict, all_good_bins, starting_completeness,
                                                        min_purity, min_completeness, threads, n_dim, annot_file,
                                                        mg_depth_file, single_contig_bins, taxon_marker_sets,
                                                        tigrfam2pfam_data, main_contig_data_dict)

all_contigs = []
for bin in all_good_bins:
    all_contigs.extend(all_good_bins[bin]['contigs'])
if len(all_contigs) != len(set(all_contigs)):
    logging.warning('WARNING: {0} duplicate contigs in bins found!'.format(len(all_contigs) - len(set(all_contigs))))

# Write bin fastas.
write_bins(all_good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity),
           bin_dir='bins')

all_cont_data_dict = {}

for contig, k_freq in main_contig_data_dict.items():
    all_cont_data_dict[contig] = {'k-mer_freqs': list(k_freq), 'depths': list(depth_dict.get(contig))}

with open("contig_data_dict.pickle", "wb") as of:
    pickle.dump(all_cont_data_dict, of, pickle.HIGHEST_PROTOCOL)

logging.info('Run finished.')
# sys.stdout.close()
