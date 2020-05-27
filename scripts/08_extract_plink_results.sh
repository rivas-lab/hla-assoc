#!/bin/bash

# Get allelotypes that are >= 0.1% frequency
output="../output/hla_plink_output.tsv"

rm $output
touch $output

echo -e "GBE_ID\tCHROM\tPOS\tALL_ID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP" > $output
cut -f1 ../hla_genotype_data/haps_above_freq_thresh.tsv | tail -n +2 > tmp

# Get phe files that are >= 500
echo "Extracting PLINK summarystatistics for HLA alleles above frequency 0.1% and phenotypes with above 500 cases..."
phe_files=$(cut -f15 ../output/hla_phenotype_info_raw.tsv | tail -n +2)
for file in $phe_files; do
    phe=$(basename $file .phe)
    echo $phe
    sumstat=/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/white_british/ukb24983_v2_hg19.${phe}.array-combined.glm.logistic.hybrid.gz;
    if [ -f $sumstat ]; then
        grep -Fwf tmp <(zcat $sumstat | awk -F'\t' -v phe=$phe '{print phe"\t"$0}' | cut -f1-14) >> $output;
    fi
done

# Adjusst p-values based on BY method per phenotype
echo "Performing Benjamini-Yekutieli adjustment across all phenotypes..."
python adjust_p.py ../output/hla_plink_output.tsv "P" "BY_ADJ_P"

rm tmp
# Extract significant associations post-adjustment
echo "Extracting BY-significant associations with FDR 0.05..."
awk -F'\t' '{if ((NR == 1) || (($15 != NA) && ($15 <= 5e-2))) {print}}' ../output/hla_plink_output_adj.tsv | sort -k15 -g > ../output/hla_plink_sig_output_adj.tsv
