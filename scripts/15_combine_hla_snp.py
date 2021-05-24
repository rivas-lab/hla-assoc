#!/bin/python
import pandas as pd

hla = pd.read_table('../output/hla_bma_sig_output.tsv')[['allelotype','GBE_ID']]
hla_map = pd.read_csv('../html/HLA_alleles/ukb_to_asterisk_names.csv')
hla = hla.merge(hla_map, left_on='allelotype', right_on='literature_names')
hla = hla[['ukb_names', 'GBE_ID']]
hla.columns = ['allelotype', 'GBE_ID']
snp = pd.read_table('../output/hla_noncoding_snps_sig_output.tsv')[['ALL_ID','GBE_ID']]
snp.columns = ['allelotype', 'GBE_ID']
merged = pd.concat([hla, snp]).sort_values('GBE_ID')
merged.to_csv('../output/hla_combined_bma_snp.tsv', sep='\t', index=False)
