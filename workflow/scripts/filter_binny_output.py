#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Fri Apr  3 08:22:12 2020

@author: oskar.hickl
'''

import sys
from pathlib import Path

input_assembly = Path(sys.argv[1])
input_binny_out = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
minimum_completeness = int(sys.argv[4])
minimum_purity = int(sys.argv[5])


def load_fasta(fasta):
    """
    Parameters
    ----------
    fasta : PosixPath/STR
        Path to input fasta.
    Returns
    -------
    Dict with contig headers as keys and sequences as single string.
    """
    sequence_dict = {}
    with open(fasta) as f:
        current_header = None
        for line in f:
            # Get header
            if line.startswith('>'):
                # Transform sub lists of sequences from previous header to one list.
                if current_header:
                    sequence_dict[current_header] = ''.join(sequence_dict[current_header])
                line = line.strip('\n ')
                line = line.strip('>')
                sequence_dict[line] = []
                current_header = line
            # Get sequences
            else:
                line = line.strip('\n ')
                sequence_dict[current_header].append(line)
        # Transform sub lists of sequences from last header to one list.
        sequence_dict[current_header] = ''.join(sequence_dict[current_header])
    return sequence_dict


def filter_binny_out(assembly_dict, binny_out, binny_dir, min_c, min_p):
    """
    Parameters
    ----------
    assembly_dict : PosixPath/STR
        Path to assembly dict with contig names as keys and contig sequences as single string.
    binny_out : PosixPath/STR
        Path to Binny contigs2bin.tsv file.
    binny_dir : PosixPath/STR
        Binny output directory path.
    min_c : INT
        Minimum estimated completeness.
    min_p : INT
        Minimum estimated purity.
    Returns
    -------
    Filter binny output, selecting only bins that pass thresholds and wirte bin fastas as well as a filtered contig2bin file.
    """
    bin_data_dict = {}
    with open(binny_out) as f:
        # Skip header
        header = next(f)
        # For each line get contig and bin it was put into.
        for line in f:
            line = line.strip('\n \t')
            contig = line.split('\t')[0]
            binny_bin = line.split('\t')[1]
            # Skip bins designated as noise, contaminated or empty.
            if not binny_bin.startswith(('N', 'B', 'E')):
                # Get bin completeness and purity estimate.
                completeness = int(binny_bin.split('_')[0].strip('C'))
                purity = int(binny_bin.split('_')[1].strip('P'))
                # Add contig/bin to dict, if it passes completeness and purity thresholds
                if completeness >= min_c and purity >= min_p:
                    if not bin_data_dict.get(binny_bin):
                        bin_data_dict[binny_bin] = [contig]
                    else:
                        bin_data_dict.get(binny_bin).append(contig)
    # Create bin folder, if it doesnt exist.
    bin_dir = binny_dir / 'bins'
    bin_dir.mkdir(parents=True, exist_ok=True)
    # Write bin fasta for each contig.
    for bin in bin_data_dict:
        bin_file_name = 'binny_' + bin + '.fasta'
        bin_out_path = bin_dir / bin_file_name
        with open(bin_out_path, 'w') as out_file:
            for contig in bin_data_dict.get(bin):
                out_file.write('>' + contig + '\n' + assembly_dict.get(contig) + '\n')
    print('Written binny bin fastas.')
    # Write filtered contig2bin.tsv.
    c2b_filt_path = binny_dir / 'contigs2bin_filtered.tsv'
    with open(c2b_filt_path, 'w') as out_file:
        for bin in bin_data_dict:
            for contig in bin_data_dict.get(bin):
                out_file.write(contig + '\t' + bin + '\n')
    print('Written filtered contig to bin file.')


input_assembly_dict = load_fasta(input_assembly)
filter_binny_out(input_assembly_dict, input_binny_out, output_dir, minimum_completeness, minimum_purity)

