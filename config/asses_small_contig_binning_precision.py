import pandas as pd
import glob


def load_gsa_mappings(gsa_mappings, load_gsam=False):
    mappings_dict = {}
    with open(gsa_mappings, 'r') as mappings:
        current_sample = None
        for line in mappings:
            line = line.strip('\n ')
            if line.startswith('@SampleID'):
                line_sample_data = line.split('_')[1:]
                current_data_set = '_'.join([line_sample_data[0], line_sample_data[2]])
                if current_data_set not in mappings_dict:
                    mappings_dict[current_data_set] = {}
                current_sample = '_'.join(line_sample_data[-2:])
                mappings_dict[current_data_set][current_sample] = {}
            elif line.strip() and not line.startswith('@') and not line.startswith('#'):
                line_data = line.split('\t')
                line_contig = line_data[0]
                line_contig_bin = line_data[1]
                # For binner output
                if '/' in line_contig_bin:
                    line_contig_bin = line_contig_bin.split('/')[-1].replace('.fasta', '')
                try:
                    line_contig_length = line_data[3]
                except IndexError:
                    line_contig_length = line_data[2]
                if load_gsam:
                    mappings_dict[current_data_set][current_sample][line_contig] = {'bin': line_contig_bin,
                                                                                    'length': int(line_contig_length)}
                else:
                    if line_contig_bin not in mappings_dict[current_data_set][current_sample]:
                        mappings_dict[current_data_set][current_sample][line_contig_bin] = {line_contig: int(line_contig_length)}
                    else:
                        mappings_dict[current_data_set][current_sample][line_contig_bin][line_contig] = int(line_contig_length)
    return mappings_dict


def load_binner_bin_metrics(binner_bin_metrics):
    bin_dict = {}
    with open(binner_bin_metrics, 'r') as bin_metrics:
        header = next(bin_metrics).strip('\n ').split('\t')[2:]
        for line in bin_metrics:
            line = line.strip('\n ')
            line_bin_data = line.split('\t')
            line_data_set_data = line_bin_data[0].split('_')
            line_data_set = '_'.join([line_data_set_data[1], line_data_set_data[3]])
            line_sample = '_'.join(line_data_set_data[-2:])
            line_bin = line_bin_data[1]
            if '/' in line_bin:
                line_bin = line_bin.split('/')[-1].replace('.fasta', '')
            bin_dict[line_bin] = {stat: (val if stat == 'Most abundant genome' else round(float(val), 4) if '.' in val
                                         else int(val))
                                  for stat, val in zip(header, line_bin_data[2:])}
            bin_dict[line_bin]['Data set'] = line_data_set
            bin_dict[line_bin]['Sample ID'] = line_sample
    return bin_dict


amber_run = 'cami_all'

gsa_mapping = '/Users/oskar.hickl/binny_bench/amber_data/{0}/gsa_mappings.tsv'.format(amber_run)
gsam_contig_data = load_gsa_mappings(gsa_mapping, load_gsam=True)

binny_version = 'Binny'
binner_mapping = '/Users/oskar.hickl/binny_bench/amber_data/{0}/{1}_bins_amber.tsv'.format(amber_run, binny_version)
bin_stats = '/Users/oskar.hickl/binny_bench/amber_data/{0}/genome/{1}/metrics_per_bin.tsv'.format(amber_run, binny_version)
binny_contig_data = load_gsa_mappings(binner_mapping)
bin_stat_data = load_binner_bin_metrics(bin_stats)

# max_length = 1000
for max_length in [1000]:  # [1000000, 3000, 2000, 1000, 500]
    df_rows = []
    # small_contigs = 0
    # correctly_assigned_contigs = 0
    for bin, bin_data in bin_stat_data.items():
        small_contigs = 0
        correctly_assigned_contigs = 0

        size_small_contigs = 0
        size_correctly_assigned_contigs = 0

        bin_size = bin_data['Bin size (bp)']
        reference_size = bin_data['True size of most abundant genome (bp)']

        bin_data_set = bin_data['Data set']
        bin_data_sample = bin_data['Sample ID']

        bin_contigs = [contig for contig in binny_contig_data[bin_data_set][bin_data_sample][bin].keys()
                       if gsam_contig_data[bin_data_set][bin_data_sample][contig]['length'] + 1 < max_length]
        if bin_contigs:
            # print(max_length, len(bin_contigs))
            for contig in bin_contigs:
                size_contig = gsam_contig_data[bin_data_set][bin_data_sample][contig]['length'] + 1
                small_contigs += 1
                size_small_contigs += size_contig
                if bin_stat_data[bin]['Most abundant genome'] == gsam_contig_data[bin_data_set][bin_data_sample][contig]['bin']:
                    correctly_assigned_contigs += 1
                    size_correctly_assigned_contigs += size_contig
            bin_row = [bin_data_set, bin_data_sample, bin, bin_size, reference_size, small_contigs, size_small_contigs,
                       round(size_small_contigs / bin_size, 4), correctly_assigned_contigs, size_correctly_assigned_contigs,
                       round(size_correctly_assigned_contigs / bin_size, 4),
                       round(correctly_assigned_contigs / small_contigs, 4)]
            df_rows.append(bin_row)

    # print(small_contigs, correctly_assigned_contigs, round(correctly_assigned_contigs / small_contigs * 100, 1))

    short_contigs_precison_df = pd.DataFrame(df_rows, columns=['Data set', 'Sample', 'Bin', 'Bin size[bp]',
                                                               'Reference size[bp]', 'Small contigs',
                                                               'Size small contigs[bp]', 'Fraction small contigs[bp]',
                                                               'tp small contigs', 'Size tp small contigs[bp]',
                                                               'Fraction tp small contigs[bp]', 'Precision'])
    short_contigs_precison_df.to_csv('/Users/oskar.hickl/binny_bench/amber_data/{0}/'
                                     '{1}_cami_contigs_smaller_{2}_precison.tsv'.format(amber_run, binny_version, max_length),
                                     sep='\t', index=False)

    # Get stats for bins, where small contigs make up more than X % of bin size in bp.
    for perc in ['025']:  # ['00', '01', '05', '10', '25']
        short_contigs_precison_df.query('`Fraction small contigs[bp]` > 0.{0}'
                                        .format(perc)).describe().to_csv('/Users/oskar.hickl/binny_bench/amber_data/{0}/'
                                                                 '{1}_cami_contigs_smaller_{2}_precision_more_than_{3}_perc_share'
                                                                 '_summary.tsv'.format(amber_run, binny_version, max_length, perc),
                                                                         sep='\t')

# max_length = 500
# n_bins = 0
# for bin, bin_data in bin_stat_data.items():
#     bin_data_set = bin_data['Data set']
#     bin_data_sample = bin_data['Sample ID']
#     # for contig in binny_contig_data[bin_data_set][bin_data_sample][bin].keys():
#     #     if gsam_contig_data[bin_data_set][bin_data_sample][contig]['length'] < max_length:
#     #         print(contig, gsam_contig_data[bin_data_set][bin_data_sample][contig]['length'])
#     bin_contigs = [contig for contig in binny_contig_data[bin_data_set][bin_data_sample][bin].keys()
#                    if gsam_contig_data[bin_data_set][bin_data_sample][contig]['length'] < max_length]
#     # for contig in bin_contigs:
#     #     print(contig, gsam_contig_data[bin_data_set][bin_data_sample][contig]['length'])
#     if bin_contigs:
#         n_bins += 1
# print(n_bins)
