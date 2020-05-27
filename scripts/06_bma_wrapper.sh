#/bin/bash

python adjust_p.py ../output/hla_additive_output.tsv "P" "BY_ADJ_P"

# Extract significant associations post-adjustment
echo "Extracting BY-significant associations with FDR 0.05..."
awk -F'\t' '{if ((NR == 1) || ($8 <= 5e-2)) {print}}' ../output/hla_additive_output_adj.tsv | sort -k8 -g > ../output/hla_additive_sig_output_adj.tsv

batch_dir="../output/bma"

phes=$(cut -f1 ../output/hla_additive_sig_output_adj.tsv | sort | uniq -c | awk '{if ($1 > 1) {print $2}}')

# for every phenotype that has more than 2 significant associations, run BMA
for pheno in $phes; do
    phe_path=$(awk -F'\t' -v pheno=$pheno '{if ($1 == pheno) {print}}' ../output/hla_phenotype_info_raw.tsv | cut -f15)
    echo $phe_path
    sbatch -J $pheno -o $batch_dir/hla_bma_${pheno}.out  -t 06:00:00 -p normal,mrivas,owners -N 1 --mem=16000 --wrap="Rscript 06_bma.R $phe_path 20 BIN"
done
