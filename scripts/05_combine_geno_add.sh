#!/bin/bash

output_dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/hla_manuscript/output/sumstats"

rm ../output/hla_genotype_output.tsv
rm ../output/hla_additive_output.tsv

echo -e "GBE_ID\tALL_ID\tBETA\tSE\tZ_STAT\tP\tAIC" > ../output/hla_genotype_output.tsv
echo -e "GBE_ID\tALL_ID\tBETA\tSE\tZ_STAT\tP\tAIC" > ../output/hla_additive_output.tsv

cut -f1 ../hla_genotype_data/haps_above_freq_thresh.tsv | tail -n +2 > tmp

for file in $(ls $output_dir/hla_genotype_*tsv); do
    phe=$(basename $file .tsv | cut -d'_' -f3);
    echo $phe;
    grep -Fwf tmp <(awk -F'\t' -v phe=$phe '{if (NR > 1) {print phe"\t"$0}}' $file | egrep -v "age|sex|Array|PC|CNV|Intercept") >> ../output/hla_genotype_output.tsv
done

for file in $(ls $output_dir/hla_additive_*tsv); do
    phe=$(basename $file .tsv | cut -d'_' -f3);
    echo $phe;
    grep -Fwf tmp <(awk -F'\t' -v phe=$phe '{if (NR > 1) {print phe"\t"$0}}' $file | egrep -v "age|sex|Array|PC|CNV|Intercept") >> ../output/hla_additive_output.tsv
done

rm tmp
sed 's/.as.numeric(gcounts)//g' ../output/hla_additive_output.tsv > tmp && mv tmp ../output/hla_additive_output.tsv
