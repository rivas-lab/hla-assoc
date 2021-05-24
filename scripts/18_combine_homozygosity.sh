#!/bin/bash

output_dir="/oak/stanford/groups/mrivas/users/guhan/repos/hla-assoc/output/homozygosity"

rm ../output/hla_homozygosity_output.tsv

echo -e "GBE_ID\tALL_ID\tBETA\tSE\tZ_STAT\tP\tAIC\tLOCUS" > ../output/hla_homozygosity_output.tsv

for file in $(ls $output_dir/hla_homozygosity_*tsv); do
    phe=$(basename $file .tsv | cut -d'_' -f3);
    echo $phe;
    awk -F'\t' -v phe=$phe '{if (NR > 1) {print phe"\t"$0}}' $file | grep gcounts >> ../output/hla_homozygosity_output.tsv
done

sed 's/dosage.gcounts/homozygote_count/g' ../output/hla_homozygosity_output.tsv >tmp && mv tmp ../output/hla_homozygosity_output.tsv

python adjust_p.py ../output/hla_homozygosity_output.tsv "P" "BY_ADJ_P"

grep -Fwf <(cut -f4 ../output/hla_complete.tsv | sort -u) ../output/hla_homozygosity_output_adj.tsv > tmp && mv tmp ../output/hla_homozygosity_output_adj.tsv
