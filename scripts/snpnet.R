#!/bin/R

args <- commandArgs(TRUE)
phenotype <- args[1]

#install.packages('/oak/stanford/groups/mrivas/users/ruilin/softwares/myglmnet', repo=NULL,type='source')
#install.packages('/oak/stanford/groups/mrivas/users/ruilin/softwares/plink-ng/2.0/pgenlibr', repo=NULL,type='source')
library(snpnet)

phenotype.file <- "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20201002.phe"
genotype.pfile = "/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3"

family <- "binomial"
results.dir <- paste0("../output/snpnet/", phenotype)
dir.create(results.dir)

if(family != "cox"){
  covariates <- c("age", "sex", paste0("PC", 1:10))
  status = NULL
} else {
  covariates <- c( "sex", paste0("PC", 1:10))
  status = paste0("coxnet_status_f.", phenotype, ".0.0")
  phenotype = paste0("coxnet_y_f.", phenotype, ".0.0")
  phenotype.file = file.path(scratch_dir, "master.cox.20200413.allinds.phe.zst")
}


configs <- list(
  results.dir = results.dir,
  nCores = 16,
  num.snps.batch = 2000,
  nlams.init = 20,
  nlams.delta = 10,
  save = TRUE,
  niter = 100,
  glmnet.thresh = 10^(-7),
  prevIter = 0,
  use.glmnetPlus = TRUE,
  early.stopping = TRUE,
  verbose = TRUE
)

start = Sys.time()
fit_snpnet <- snpnet(
  genotype.pfile = genotype.pfile,
  phenotype.file = phenotype.file,
  phenotype = phenotype,
  covariates = covariates,
  family = family,
  split.col = "split",
  status.col = status,
  mem = 32000,
  configs = configs
)
duration = Sys.time() - start

best_ind = which.max(fit_snpnet$metric.val)
snps_selected_inorder = character(0)

for(i in 1:best_ind){
  snps_selected_inorder = union(snps_selected_inorder, names(which(fit_snpnet$beta[[i]] != 0)))
}
lapply(snps_selected_inorder, write, paste0(results.dir, "/", phenotype, "_snpnet.txt"), append=TRUE)
