#!/bin/python

import pandas as pd

# Phenotypes with more than 500 WB individuals for which analyses have been run
gbe = pd.read_table('../output/hla_phenotype_info_raw.tsv')

# Deduplicate phenotypes. First, remove exact duplicates
gbe['CLEAN_NAME'] = gbe['GBE_NAME'].str.replace('TTE_', '').str.replace('AD_', '')
gbe = gbe.sort_values('CLEAN_NAME')
duplicated = list(gbe[gbe['CLEAN_NAME'].duplicated()]['CLEAN_NAME'])
print("Exact duplicates:")
print(duplicated)
gbe = gbe.sort_values('N_GBE', ascending=False)
idx = gbe.groupby(['CLEAN_NAME'])['N_GBE'].transform(max) == gbe['N_GBE']
clean_gbe = gbe[idx].sort_values('CLEAN_NAME')

# These are the synonyms I could find
to_remove = ['HC171','HC132','HC415','HC222','HC281','HC440','HC371','HC188','HC442','HC322','HC316','HC273','HC6','HC209','HC164','HC991','HC1574','HC299','HC226','HC198','HC287','HC311','HC45','HC123','HC156','HC285','HC370','HC206','HC54','HC146','HC203','HC1061','HC302','HC149','HC224','HC397','HC352','HC165','HC43','HC476','HC55','HC643','HC1213']
print("Fuzzy + manual curation duplicates:" + str(len(set(to_remove))))

clean_gbe = clean_gbe[~clean_gbe['#GBE_ID'].isin(to_remove)]
clean_gbe.fillna("NA").drop(columns=['CLEAN_NAME']).to_csv('../output/hla_phenotype_info.tsv', sep='\t', index=False)

complete = pd.read_table('../output/hla_complete_raw.tsv')
complete = complete[complete['GBE_ID'].isin(clean_gbe['#GBE_ID'])]

print(len(complete['GBE_ID'].unique()))
print(complete.isna().sum())
complete.to_csv('../output/hla_complete.tsv', sep='\t', index=False)

complete = complete[complete['BMA_posterior_prob'] >= 80]
complete.to_csv('../output/hla_bma_sig_output.tsv', sep='\t', index=False)
