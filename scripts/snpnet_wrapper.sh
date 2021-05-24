#/bin/bash

batch_dir="../output/snpnet"

phes=$(cut -f1 ../output/hla_additive_sig_output_adj.tsv | sort | uniq -c | awk '{if ($1 > 1) {print $2}}')

# for every phenotype that has more than 2 significant associations, run BMA
for pheno in $phes; do
    echo $pheno
    sbatch -J $pheno -o $batch_dir/hla_snpnet_${pheno}.out  -t 01:00:00 -p normal,mrivas,owners -N 1 --mem=16000 --wrap="ml load R/3.6.1 && ml msc && Rscript snpnet.R $pheno"
done
