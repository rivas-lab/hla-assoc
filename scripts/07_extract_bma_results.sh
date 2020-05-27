#!/bin/bash

output_dir="/oak/stanford/groups/mrivas/users/guhan/sandbox/hla_manuscript/output/bma"

rm ../output/hla_bma_output.tsv

echo -e "ALL_ID\tposterior_mean\tposterior_sd\tposterior_prob\tGBE_ID" > ../output/hla_bma_sig_output.tsv

for file in $(ls $output_dir | grep tsv); do
    egrep -v "age|sex|Array|PC|CNV" $output_dir/$file | awk -F'\t' '{if (NR > 1) {print}}' >> ../output/hla_bma_output.tsv;
done
