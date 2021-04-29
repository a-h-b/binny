
import sys
from pathlib import Path
import shutil

proj_path = Path('/mnt/lscratch/users/ohickl/binning/binny_eval')

binner = 'maxbin'

runs = list(proj_path.glob('{0}_out/2017.12.04_18.45.54_sample_?'.format(binner)))
runs += list(proj_path.glob('{0}_out/2017.12.04_18.45.54_sample_??'.format(binner)))

good_bins_dict = {str(sample).split('/')[-1]: [] for sample in runs}

for run in runs:
    checkm_out = run / 'checkm/output.txt'
    sample = str(run).split('/')[-1]
    if checkm_out.exists():
        with open(checkm_out, 'r') as co:
            for line in co:
                if not line.startswith(('-', '  Bin Id')):
                    line = line.split()
                    if int(float(line[12])) >= 70 and int(float(line[13])) <= 10:
                        # print(line)
                        good_bins_dict[sample].append(line[0])
# print(len(good_bins_dict['H_S001']))

for k, v in good_bins_dict.items():
    out_folder = Path('/mnt/lscratch/users/ohickl/binning/binny_eval/{0}_out/{1}/bins_checkm'.format(binner, k))
    out_folder.mkdir(exist_ok=True)
    for bin in v:
        bin = bin + '.fasta'
        org_file = Path('/mnt/lscratch/users/ohickl/binning/binny_eval/{0}_out/{1}/bins/{2}'.format(binner, k, bin))
        dest_file = out_folder / bin
        shutil.copy(org_file, dest_file)

