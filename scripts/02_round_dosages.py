#!/bin/python

import numpy as np
import time
import pandas as pd


# rounds value to nearest integer if within threshold; otherwise returns -1
def round_dosages(val, thresh):
    if val < thresh:
        return 0
    elif 1 - thresh < val and val < 1 + thresh:
        return 1
    elif 2 - thresh < val:
        return 2
    else:
        return -1

t0 = time.time()

dosage_file_name = '/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v2.txt'
fam_file_name = '/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.fam'

print('Reading in dosage file: time {}'.format(time.time() - t0))
dosage_values = np.loadtxt(dosage_file_name, skiprows = 1)

print('Reading in fam file: time {}'.format(time.time() - t0))
fam_info = np.loadtxt(fam_file_name, dtype = 'string_', usecols = range(5))

print('Rounding: time {}'.format(time.time() - t0))
thresh = 0.1
v_round_dosages = np.vectorize(round_dosages)
rounded_dosages = v_round_dosages(dosage_values, thresh)

print('Writing to file: time {}'.format(time.time() - t0))
all_haps = list(pd.read_csv('../hla_genotype_data/all_haps.tsv', sep='\t')['ID'])
df = pd.DataFrame(data=rounded_dosages, columns=all_haps)
df.to_csv("/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v2_rounded.txt",  sep='\t', index=False)
