#!/bin/R

library("dplyr")
library("stats")
library("data.table")

args <- commandArgs(TRUE)
phe_path <- args[1]
options(warn=1)
phe_name <- strsplit(phe_path, "/")[[1]]
phe_name <- strsplit(phe_name[length(phe_name)], "[.]")[[1]][1]

all_ids = as.data.frame(read.table("../hla_genotype_data/haps_above_freq_thresh.tsv", header=TRUE))

print("Reading in dosage file...")
dosage <- readRDS("../hla_genotype_data/ukb_hla_v3_24983_rounded_wb.rds")

print("Subsetting to analysis alleles...")
dosage <- dosage[, names(dosage) %in% all_ids$ID]

print("Reading in covariates file...")
covars <- readRDS("../hla_genotype_data/ukb_hla_v3_covar_24983_wb.rds")

print("Reading phenotype:")
print(phe_name)
print("from file:")
print(phe_path)
phe <- as.data.frame(read.table(phe_path, header=FALSE, row.names=1))
subset_ids <- intersect(rownames(covars), rownames(phe))

print("Subsetting to those rows that have covariates...")
phe <- phe[subset_ids,]
covars <- covars[subset_ids,]
dosage <- dosage[subset_ids,]

print("Adding phenotype to covariate matrix, filtering out missing phenotypes...")
covars["status"] = as.numeric(phe$V3)
covars <- subset(covars, status!=-9)

print("Replacing binary values...")
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

print("Subsetting dosage file to those rows in covariates file...")
if (dim(covars)[1] != dim(dosage)[1]) {
    dosage <- dosage[rownames(covars),]
}

input = cbind(dosage, covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "N_CNV", "LEN_CNV", "status")])

print("Starting GWAS...")
results_gen <- list()
results_add <- list()

for (allele in colnames(dosage)) {
    print(allele)
    # Nonadditive
    input["gcounts"] <- as.numeric(input[,allele])
    fit  <- glm(status ~ as.factor(gcounts) + scale(age) + sex + Array + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +  scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(N_CNV) + scale(LEN_CNV),
           family="binomial", data=input)
    results_df <- as.data.frame(summary.glm(fit)$coefficients)
    results_df$aic <- fit$aic
    #rownames(results_df) <- paste(allele, "_", rownames(results_df), sep="")
    results_gen[[allele]] <- results_df
    # Additive
    fit <- glm(status ~ as.numeric(gcounts) + scale(age) + sex + Array + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +  scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(N_CNV) + scale(LEN_CNV),
           family="binomial", data=input)
    results_df <- as.data.frame(summary.glm(fit)$coefficients)
    results_df$aic <- fit$aic
    #rownames(results_df) <- paste(allele, "_", rownames(results_df), sep="")
    results_add[[allele]] <- results_df
}

print("Merging tables...")
gen <- do.call(rbind, results_gen)
gen <- setDT(gen, keep.rownames="ALL_ID")
add <- do.call(rbind, results_add)
add <- setDT(add, keep.rownames="ALL_ID")

print("Writing to file...")
write.table(add, paste("../output/sumstats/hla_additive_", phe_name, ".tsv", sep=""), sep='\t', quote=FALSE, row.names=FALSE)
write.table(gen, paste("../output/sumstats/hla_genotype_", phe_name, ".tsv", sep=""), sep='\t', quote=FALSE, row.names=FALSE)
