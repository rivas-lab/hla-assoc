#!/bin/R
library(survival)
library(dplyr)
library(pgenlibr)
library(data.table)

genotype.pfile="/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3"
phenotype.file=as.data.frame(read.table("/oak/stanford/groups/mrivas/users/justesen/projects/disease_progression/phe/snpnet_cox.f.40007.0.0_new.phe", header=TRUE, row.names=1))
status = "coxnet_status_f.40007.0.0"
phenotype = "coxnet_y_f.40007.0.0"
covariates <- c("Array", "sex", paste0("PC", 1:10))
family="cox"

wb <- as.data.frame(read.table("/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe", header=FALSE))
colnames(wb) <- c("FID", "IID")
wb_ids <- as.numeric(unlist(wb %>% select("IID")))

psam <- as.data.frame(read.table("/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.psam", header=FALSE))
colnames(psam) <- c("FID", "IID", "status")
psam_ids <- as.numeric(unlist(psam %>% select("IID")))

pvar <- as.data.frame(read.table("/oak/stanford/groups/mrivas/ukbb24983/hla/pgen/ukb_hla_v3.pvar", header=FALSE))
colnames(pvar) <- c("CHROM", "POS", "ID", "REF", "ALT")
pvar_vars <- as.vector(unlist(pvar %>% select("ID")))

pgen = pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), sample_subset=match(wb_ids, psam_ids))
datamatrix = as.data.frame(pgenlibr::ReadList(pgen, variant_subset=1:362, meanimpute=T))
colnames(datamatrix) <- pvar_vars
row.names(datamatrix) <- wb_ids

input <- merge(datamatrix, phenotype.file, by=0)

form = as.formula(paste("Surv(", phenotype, ",", status, ") ~ ", paste(c(covariates, pvar_vars), collapse = " + ")))
coxmod = coxph(form, data=input, control=coxph.control(iter.max = 50))

cox_sum <- summary(coxmod)$coefficients
print(cox_sum)
write.table(cox_sum, file="../output/coxph.tsv", quote=F, sep="\t")
