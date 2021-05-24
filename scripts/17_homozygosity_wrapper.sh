#!/bin/bash

batch_dir="../output/homozygosity"

phes=$(cut -f15 ../output/hla_phenotype_info_raw.tsv | tail -n +2)

for phe in $phes; do
    echo $(basename $phe .phe);
    sbatch -J $(basename $phe .phe) -o $batch_dir/$(basename $phe .phe).na.out -t 02:00:00 -p normal,mrivas,owners -N 1 --mem=16000 --wrap="Rscript 17_homozygosity.R $phe";
done
