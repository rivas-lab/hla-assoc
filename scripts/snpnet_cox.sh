#!/bin/bash
#SBATCH --job-name=snpnet_cox
#SBATCH --output=logs/snpnet_cox.%A.out
#SBATCH --error=logs/snpnet_cox.%A.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=6000
#SBATCH --time=01:00:00
#SBATCH -p normal,owners,mrivas
#SBATCH -C "CPU_GEN:HSW|CPU_GEN:BDW|CPU_GEN:SKX"
set -beEuo pipefail
ml msc
############################################################
# Required arguments for ${snpnet_wrapper} script
############################################################
genotype_pfile="/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3"
phe_file="/oak/stanford/groups/mrivas/users/justesen/projects/disease_progression/phe/snpnet_cox.f.40007.0.0_new.phe"
phenotype_name="coxnet_y_f.40007.0.0" # One may use phenotype_name=$1 etc
family="cox"
results_dir="/oak/stanford/groups/mrivas/users/guhan/repos/hla-assoc/output/snpnet_cox/${phenotype_name}"

snpnet_dir="/oak/stanford/groups/mrivas/users/justesen/projects/snpnet"
snpnet_wrapper="/oak/stanford/groups/mrivas/users/justesen/projects/snpnet/helpers/snpnet_wrapper.sh"
############################################################
# Additional optional arguments for ${snpnet_wrapper} script
############################################################
covariates="sex,Array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
split_col="split"
status_col="coxnet_status_f.40007.0.0"

############################################################
# Configure other parameters
############################################################
cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
#ml load snpnet
# Two variables (${snpnet_dir} and ${snpnet_wrapper}) should be already configured by Sherlock module
# https://github.com/rivas-lab/sherlock-modules/tree/master/snpnet
# Or, you may use the latest versions
#  snpnet_dir="$OAK/users/$USER/repos/rivas-lab/snpnet"
#  snpnet_wrapper="$OAK/users/$USER/repos/rivas-lab/PRS/helper/snpnet_wrapper.sh"

############################################################
# Run ${snpnet_wrapper} script
############################################################

echo "[$0 $(date +%Y%m%d-%H%M%S)] [start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; phenotype = ${phenotype_name}" >&2

bash ${snpnet_wrapper} \
--snpnet_dir ${snpnet_dir} \
--nCores ${cores} --memory ${mem} \
--covariates ${covariates} \
--split_col ${split_col} \
--status_col ${status_col} \
--verbose \
${genotype_pfile} \
${phe_file} \
${phenotype_name} \
${family} \
${results_dir}

# --no_save

echo "[$0 $(date +%Y%m%d-%H%M%S)] [end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}; phenotype = ${phenotype_name}" >&2

############################################################
# Another example (w/ the sample data in snpnet package)
############################################################
#genotype_pfile="$OAK/users/$USER/repos/rivas-lab/snpnet/inst/extdata/sample"
#phe_file="$OAK/users/$USER/repos/rivas-lab/snpnet/inst/extdata/sample.phe"
#phenotype_name="QPHE"
#family="gaussian"
#results_dir="$OAK/users/$USER/repos/rivas-lab/PRS/notebook/20191021_snpnet/private_out/20/${phenotype_name}"
#results_dir="/scratch/users/ytanigaw/snpnet.demo/${phenotype_name}"
