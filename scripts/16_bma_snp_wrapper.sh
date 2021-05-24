#/bin/bash

batch_dir="../output/bma_snp"

phes=$(cut -f2 ../output/hla_combined_bma_snp.tsv | tail -n +2 | sort | uniq -c | awk '{if ($1 > 1) {print $2}}')

# for every phenotype that has more than 2 significant HLA/SNP associations, run BMA (77)
for pheno in $phes; do
    phe_path=$(awk -F'\t' -v pheno=$pheno '{if ($1 == pheno) {print}}' ../output/hla_phenotype_info_raw.tsv | cut -f15)
    echo $phe_path
    sbatch -J $pheno -o $batch_dir/hla_bma_snp_${pheno}.out  -t 06:00:00 -p normal,mrivas,owners -N 1 --mem=16000 --wrap="Rscript 16_bma_snp.R $phe_path 20 BIN"
done
