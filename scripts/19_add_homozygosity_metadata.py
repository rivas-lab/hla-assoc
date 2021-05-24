#!/bin/python

import pandas as pd
import numpy as np

phen_info = pd.read_table('../output/hla_phenotype_info_raw.tsv', dtype=object)[['#GBE_ID', 'GBE_NAME']]
phen_info.columns = ['GBE_ID', 'GBE_NAME']

hz = pd.read_table('../output/hla_homozygosity_output_adj.tsv')

err = 1.96 * hz['SE'].astype(float)
hz['u95_BETA'] = hz['BETA'].astype(float) + err
hz['l95_BETA'] = hz['BETA'].astype(float) - err
hz['OR'] = np.exp(hz['BETA'].astype(float))
hz['CI'] = hz['OR'].round(2).astype(str) + " [" + np.exp(hz['l95_BETA'].astype(float)).round(2).astype(str) + ", " + np.exp(hz['u95_BETA'].astype(float)).round(2).astype(str) + "]"
print(len(hz))
hz = hz.merge(phen_info).sort_values('BY_ADJ_P').drop(columns='ALL_ID')
print(len(hz))
hz.to_csv('../output/hla_homozygosity_output_adj_complete.tsv', sep='\t', index=False)
