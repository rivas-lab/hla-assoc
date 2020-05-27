#!/bin/python

import pandas as pd
import numpy as np
from collections import defaultdict
from adjust_p import BY_adjust_p

all_plink = pd.read_table('../output/hla_plink_output_adj.tsv', dtype=object).drop(columns=['REF','ALT','A1','FIRTH?','TEST','CHROM','POS'])
all_plink.columns = ['GBE_ID', 'ALL_ID', 'plink_OBS_CT', 'plink_OR', 'plink_LOG(OR)_SE', 'plink_Z_STAT', 'plink_P', 'plink_BY_ADJ_P']

all_bma = pd.read_table('../output/hla_bma_output.tsv', dtype=object)
all_bma.columns = ['ALL_ID', 'BMA_posterior_mean', 'BMA_posterior_sd', 'BMA_posterior_prob', 'GBE_ID']

all_geno = pd.read_table('../output/hla_genotype_output.tsv', dtype=object)
data = defaultdict(dict)
for i, row in all_geno.iterrows():
    gen = row['ALL_ID'].split('.as.factor(gcounts)')[1]
    all_id = row['ALL_ID'].replace('.as.factor(gcounts)' + gen,'')
    data[(row['GBE_ID'], all_id)]['gen_' + gen + '_BETA'] = row['BETA']
    data[(row['GBE_ID'], all_id)]['gen_' + gen + '_SE'] = row['SE']
    data[(row['GBE_ID'], all_id)]['gen_' + gen + '_Z_STAT'] = row['Z_STAT']
    data[(row['GBE_ID'], all_id)]['gen_' + gen + '_P'] = row['P']
    data[(row['GBE_ID'], all_id)]['AIC_genotype'] = row['AIC']
all_geno = pd.DataFrame(data).T.reset_index()
all_geno = all_geno.rename(columns={'level_0': 'GBE_ID', 'level_1': 'ALL_ID'})

all_add = pd.read_table('../output/hla_additive_output_adj.tsv', dtype=object)
all_add.columns = ['GBE_ID', 'ALL_ID', 'add_BETA', 'add_SE', 'add_Z_STAT', 'add_P', 'AIC_additive', 'add_BY_ADJ_P']
all_add['ALL_ID'] = all_add['ALL_ID'].str.replace('.as.numeric\(gcounts\)','')

phen_info = pd.read_table('../output/hla_phenotype_info_raw.tsv', dtype=object)[['#GBE_ID', 'N_GBE', 'GBE_NAME']]
phen_info.columns = ['GBE_ID', 'N_GBE', 'GBE_NAME']

merged = all_plink.merge(all_bma, how='outer', on=['GBE_ID','ALL_ID'])
print(merged.isna().sum())
merged = merged.merge(all_geno, how='outer', on=['GBE_ID','ALL_ID'])
print(merged.isna().sum())
merged = merged.merge(all_add, how='outer', on=['GBE_ID','ALL_ID'])
print(merged.isna().sum())
merged = merged.merge(phen_info, how='left', on=['GBE_ID'])
merged['delta_AIC'] = merged['AIC_genotype'].astype(float) - merged['AIC_additive'].astype(float)
# 698 sumstats not generated across R and PLINK both
merged = merged[merged['AIC_additive'].notna()]

merged = merged.merge(pd.read_csv('/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_to_asterisk_names.csv'), left_on='ALL_ID', right_on='ukb_names')
all_haps = pd.read_csv('../hla_genotype_data/all_haps.tsv', sep='\t')
all_haps.columns = ['ALL_ID', 'frequency']
merged = merged.merge(all_haps, how='left')

merged['allelotype'] = merged['literature_names']
print(merged.isna().sum())
# 7620 records of allelotypes only having a max of 1 allele count, no 2

merged['gen_1_P'] = merged['gen_1_P'].astype(float)
merged['gen_2_P'] = merged['gen_2_P'].astype(float)
merged = BY_adjust_p(merged, 'gen_1_P', 'gen_1_BY_ADJ_P')
merged = BY_adjust_p(merged, 'gen_2_P', 'gen_2_BY_ADJ_P')
merged['locus'] = ('HLA-' + merged['allelotype'].apply(lambda x: x.split('*')[0]))

for col in ['add','gen_1','gen_2']:
    err = 1.96 * merged['{}_SE'.format(col)].astype(float)
    merged['u95_' + col + '_BETA'] = merged[col + '_BETA'].astype(float) + err
    merged['l95_' + col + '_BETA'] = merged[col + '_BETA'].astype(float) - err
    merged[col + '_OR'] = np.exp(merged[col + '_BETA'].astype(float))
    merged[col + '_CI'] = merged[col + '_OR'].round(2).astype(str) + " [" + np.exp(merged['l95_' + col + '_BETA'].astype(float)).round(2).astype(str) + ", " + np.exp(merged['u95_' + col + '_BETA'].astype(float)).round(2).astype(str) + "]"

merged['u95_OR_BMA'] = np.exp(merged['BMA_posterior_mean'].astype(float) + (1.96 * merged['BMA_posterior_sd'].astype(float)))
merged['l95_OR_BMA'] = np.exp(merged['BMA_posterior_mean'].astype(float) - (1.96 * merged['BMA_posterior_sd'].astype(float)))
merged['BMA_CI'] = np.exp(merged['BMA_posterior_mean'].astype(float)).round(2).astype(str) + " [" + merged['l95_OR_BMA'].round(2).astype(str) + ", " + merged['u95_OR_BMA'].round(2).astype(str) + "]"

merged[['allelotype', 'locus', 'frequency', 'GBE_ID', 'GBE_NAME', 'N_GBE', 'plink_OR', 'plink_LOG(OR)_SE', 'plink_P', 'plink_BY_ADJ_P', 'plink_Z_STAT', 'plink_OBS_CT', 'BMA_posterior_mean', 'BMA_posterior_sd', 'BMA_posterior_prob', 'BMA_CI', 'gen_1_BETA', 'gen_1_SE', 'gen_1_P', 'gen_1_BY_ADJ_P', 'gen_1_Z_STAT', 'gen_1_CI', 'gen_2_BETA', 'gen_2_SE', 'gen_2_P', 'gen_2_BY_ADJ_P', 'gen_2_Z_STAT', 'gen_2_CI', 'add_BETA', 'add_SE', 'add_P', 'add_BY_ADJ_P', 'add_Z_STAT', 'add_CI', 'AIC_genotype', 'AIC_additive', 'delta_AIC']].sort_values(['allelotype', 'GBE_ID']).to_csv('../output/hla_complete_raw.tsv', sep='\t', index=False)
