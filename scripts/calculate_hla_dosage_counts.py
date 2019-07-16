import pandas as pd

hla = pd.read_table('/oak/stanford/groups/mrivas/users/jolivier/repos/hla-assoc/scripts/output/make_dosage_rds/ukb_hla_v2_wb_rounded_remove.tsv').drop('ID', axis = 1)
result = hla.apply(pd.value_counts, dropna=False).fillna("NA")

result.to_csv('hla_allele_dosage_counts.tsv', sep='\t')
