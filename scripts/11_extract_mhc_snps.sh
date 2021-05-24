#!/bin/bash

# Get allelotypes that are >= 0.1% frequency
output="../output/hla_snps_output.tsv"

rm $output
touch $output

echo -e "GBE_ID\tCHROM\tPOS\tALL_ID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP" > $output

# Get phe files that are >= 500
echo "Extracting PLINK summarystatistics for HLA alleles above frequency 0.1% and phenotypes with above 500 cases..."
phe_files=$(cut -f15 ../output/hla_phenotype_info_raw.tsv | tail -n +2)
for file in $phe_files; do
    phe=$(basename $file .phe)
    echo $phe
    sumstat=/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/white_british/ukb24983_v2_hg19.${phe}.array-combined.glm.logistic.hybrid.gz;
    if [ -f $sumstat ]; then
        grep -Fvwf tmp <(zcat $sumstat | awk -F'\t' -v phe=$phe '{print phe"\t"$0}' | cut -f1-14 | awk -F'\t' '{if (($2 == "6") && ($3 <= 33480577) && ($3 >= 28510120)) {print}}' | awk '$4 !~ /_+/ && $4 !~ /_-/') >> $output;
    fi
done

# Adjust p-values based on BY method per phenotype
awk -F'\t' '{if ((NR == 1) || (($14 != NA) && ($14 <= 5.887e-5))) {print}}' $output | sort -k14 -g > ../output/hla_snps_sig_output_adj.tsv
