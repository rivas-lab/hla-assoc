import pandas as pd

ann = pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv')
df = pd.read_table('../output/hla_snps_sig_output.tsv')
bma = pd.read_table('../output/hla_bma_sig_output.tsv')
miss = pd.read_table('../hla_genotype_data/missingness.tsv')

print(len(df))
df.insert(
    loc=0,
    column="V",
    value=df["CHROM"]
    .astype(str)
    .str.cat(df["POS"].astype(str), sep=":")
    .str.cat(df["REF"], sep=":")
    .str.cat(df["ALT"], sep=":"),
)

df = df.merge(ann, how='left', on='V')
intron = [
        "regulatory_region_variant",
        "intron_variant",
        "intergenic_variant",
        "downstream_gene_variant",
        "mature_miRNA_variant",
        "non_coding_transcript_exon_variant",
        "upstream_gene_variant",
        "NA",
        "NMD_transcript_variant",
]
df = df[df.most_severe_consequence.isin(intron)]
df = df[df.GBE_ID.isin(list(bma['GBE_ID']))]
df = df[(df.maf >= 0.001) & (df.ld_indep ==True)]
print(df.columns)
df = df.merge(miss, left_on="ALL_ID", right_on="ID", how="left")
df = df[df['f_miss'] < 0.001]
print(len(df))
bma_sig_nos = bma[['allelotype', 'GBE_ID']].groupby('GBE_ID').count().reset_index()

bma_sig_nos['num_snps'] = 20-bma_sig_nos['allelotype']
print(bma_sig_nos)


final_df = pd.DataFrame()

for i, row in bma_sig_nos.iterrows():
    gbe = row['GBE_ID'] 
    subset_df = df[df.GBE_ID == gbe]
    subset_df = subset_df.sort_values('P').head(row['num_snps'])
    final_df = pd.concat([final_df, subset_df])

final_df.to_csv('../output/hla_noncoding_snps_sig_output.tsv', sep='\t', index=False)
