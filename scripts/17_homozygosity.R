#!/bin/R

library("tidyverse")
library("stats")
library("data.table")

args <- commandArgs(TRUE)
phe_path <- args[1]
options(warn=1)
phe_name <- strsplit(phe_path, "/")[[1]]
phe_name <- strsplit(phe_name[length(phe_name)], "[.]")[[1]][1]

all_ids = as.data.frame(read.table("../hla_genotype_data/haps_above_freq_thresh.tsv", header=TRUE))
all_ids = all_ids %>% separate(ID, c("LOCUS", "ALLELE"), sep="_", remove=FALSE)

print("Reading in dosage file...")
dosage <- readRDS("../hla_genotype_data/ukb_hla_v3_24983_rounded_wb.rds")

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

print("Adding phenotype to covariate matrix, filtering out missing phenotypes...")
covars["status"] = as.numeric(phe$V3)
covars <- subset(covars, status!=-9)

print("Replacing binary values...")
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

count = 0
data <- list()

#Run one test per locus
for (i in unique(all_ids$LOCUS)) {
    count = count + 1
    print("Subsetting to analysis alleles...")
    subset_haps = all_ids[all_ids$LOCUS == i,]
    dosage_subset <- dosage[, names(dosage) %in% subset_haps$ID]
    
    dosage_subset$sum <- rowSums(dosage_subset == 2, na.rm=T)
    dosage_subset$nasums <-rowSums(is.na(dosage_subset))
    dosage_subset$missingness <- dosage_subset$nasums / dim(dosage_subset)[[1]]
    dosage_subset <- dosage_subset %>% select(sum,nasums,missingness)
    dosage_subset <- dosage_subset %>% select(sum,missingness)
    dosage_subset <- dosage_subset[subset_ids,]
    input = cbind(dosage_subset, covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "N_CNV", "LEN_CNV", "status")])

    print("Starting GWAS...")
    results_add <- list()

    input["gcounts"] <- as.numeric(input[,"sum"])
    print(colSums(is.na(input)))

    # Additive
    fit <- glm(status ~ scale(gcounts) + scale(age) + sex + Array + missingness + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) +  scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + scale(N_CNV) + scale(LEN_CNV),
       family="binomial", data=input)
    results_df <- as.data.frame(summary.glm(fit)$coefficients) 
    results_df$aic <- fit$aic
    results_df$locus <- i
    #rownames(results_df) <- paste(allele, "_", rownames(results_df), sep="")
    results_add[["dosage"]] <- results_df

    print("Merging tables...")
    add <- do.call(rbind, results_add)
    add <- setDT(add, keep.rownames="ALL_ID")

    data[[count]] = add
}

final <- do.call(rbind, data)
print("Writing to file...")
write.table(final, paste("../output/homozygosity/hla_homozygosity_", phe_name, ".tsv", sep=""), sep='\t', quote=FALSE, row.names=FALSE)
