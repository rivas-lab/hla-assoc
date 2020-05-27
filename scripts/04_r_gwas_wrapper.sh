#!/bin/bash

batch_dir="../output/sumstats"

grep -Fwf <(awk -F'\t' '{if ((NR == 1) || ($8 >= 500)) {print}}' /oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/05_gbe/phenotype_info.tsv | cut -f15 | egrep 'HC|cancer|PATH') /oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/05_gbe/phenotype_info.tsv > ../output/hla_phenotype_info_raw.tsv

phes=$(cut -f15 ../output/hla_phenotype_info_raw.tsv | tail -n +2)

for phe in $phes; do
    echo $(basename $phe .phe);
    sbatch -J $(basename $phe .phe) -o $batch_dir/$(basename $phe .phe).na.out -t 02:00:00 -p normal,mrivas,owners -N 1 --mem=16000 --wrap="Rscript 04_r_gwas.R $phe";
done
