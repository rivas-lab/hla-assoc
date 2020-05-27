#!/bin/bash

echo "Extracting all HLA alleles and frequencies..."
zstdcat /oak/stanford/groups/mrivas/ukbb24983/hla/afreq/ukb_hla_v3.white_british.afreq.zst | awk -F'\t' '{print $2"\t"$5}' > ../hla_genotype_data/all_haps.tsv
echo "Extracting all HLA alleles and frequencies above 0.1%..."
zstdcat /oak/stanford/groups/mrivas/ukbb24983/hla/afreq/ukb_hla_v3.white_british.afreq.zst | awk -F'\t' '{if ($5 >= 0.001) {print $2"\t"$5}}' > ../hla_genotype_data/haps_above_freq_thresh.tsv
