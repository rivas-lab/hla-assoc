args <- commandArgs(TRUE)
phe <- args[1]
phe_name <- strsplit(phe, "/")[[1]]
phe_name <- strsplit(phe_name[length(phe_name)], "[.]")[[1]][1]
print(paste("phe name", phe_name))
print("start")
options(warn=1)

print_gene <- F
t0 <- proc.time()

ukbb_files <- "../data/ukbb_files/"

# read in dosage and covars files
#dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_rounded_remove.rds")
dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_remove.rds")
covars <- readRDS("output/make_dosage_rds/ukb_hla_v2_covar_16698_remove.rds")

# read in phe file
phe <- as.data.frame(read.table(phe, header=FALSE, row.names=1))

# limit phe to individuals we're considering
phe <- phe[rownames(covars),]

print("done loading")

# use only columns with more than 5 nonzero entries
col_idx = colSums(dosage != 0) > 5
dosage <- dosage[,col_idx]
print(paste0("number excluded (additive):",362 - ncol(dosage)))

# choose which columns to include
covars <- covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4")]

# add column for case vs control
covars["status"] = as.numeric(phe$V3)
covars <- covars[is.na(covars$status) == FALSE,]

# change no case to 0, change case to 1 (from 1 and 2)
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

# remove any rows we ended up excluding before
if (dim(covars)[1] != dim(dosage)[1]) {
	dosage <- dosage[rownames(covars),]
}

# use rounded dosage instead
dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_rounded_remove.rds")
if (dim(covars)[1] != dim(dosage)[1]) {
    dosage <- dosage[rownames(covars),]
}



# get names of haplotypes that were found significant
sig_haps_df <- as.data.frame(read.table("../notebooks/output/compare_cutoff_pvals/BMA_haps_100.csv", header=FALSE, sep = ",", row.names=1, stringsAsFactors=FALSE))

hap_names <-  as.character(sig_haps_df[phe_name,])
hap_names <- hap_names[!hap_names == "NA"]
hap_names <- hap_names[!hap_names == ""]
print(hap_names)
print(paste0("hap names: ", hap_names))
print(paste("length hap names:",length(hap_names)))

results.inter <- vector("list", (length(hap_names)^2 - length(hap_names))/2)

inter_names = character((length(hap_names)^2 - length(hap_names))/2)
count <- 1
for (i in 1:(length(hap_names) - 1)) {
#    print(paste("hap_names[i]:",hap_names[i]))
    for (j in (i + 1):length(hap_names)) {
        inter_names[count] <- paste0(hap_names[i],"-",hap_names[j])
#        print(paste("hap_names[j]:",hap_names[j]))

        print(paste("curr pair:",hap_names[i],hap_names[j]))
        covars["gcounts1"] <- as.numeric(dosage[,hap_names[i]])
        covars["gcounts2"] <- as.numeric(dosage[,hap_names[j]])
        fit <- summary(glm(status ~ as.numeric(gcounts1) * as.numeric(gcounts2) + age + sex + Array + PC1 + PC2 + PC3 + PC4, family="binomial", data=covars))
        fit$deviance.resid <- NULL
        results.inter[[count]] <- fit
        count <- count + 1

    }
}

names(results.inter) <- inter_names

saveRDS(results.inter, paste0("output/interact/",phe_name,"_rounded_inter.rds"))

