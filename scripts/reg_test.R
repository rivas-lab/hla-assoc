args <- commandArgs(TRUE)
phe <- args[1]
outname <- args[2]
num_haps <- args[3]
phe_name <- strsplit(phe, "/")[[1]]
phe_name <- strsplit(phe_name[length(phe_name)], "[.]")[[1]][1]
print(paste("phe name", phe_name))
print("start")
options(warn=1)

test <- F
run_rounded <- F
run_dosage <- F
rounded_add <- F
rounded_factor <- F
run_BMA <- T
print_gene <- F
t0 <- proc.time()
#ukbb_files <- "$SCRATCH/ukbb_files/"
ukbb_files <- "../data/ukbb_files/"
if (test) {
    dosage <- readRDS("output/make_dosage_rds/test_100_remove.rds")
    covars <- readRDS("output/make_dosage_rds/test_100_covar_16698_remove.rds")
#    outname <- "test_out_dosage"
} else {
    dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_remove.rds")
    covars <- readRDS("output/make_dosage_rds/ukb_hla_v2_covar_16698_remove.rds")
#    outname <- "out_dosage"
}

phe <- as.data.frame(read.table(phe, header=FALSE, row.names=1))

# limit phe to individuals we're considering
phe <- phe[rownames(covars),]

print("done loading")

# use only columns with more than 5 nonzero entries
col_idx = colSums(dosage != 0) > 5
dosage <- dosage[,col_idx]
#sink("haps_used.txt")
#colnames(dosage)
#sink()
#print(paste0("column names: ", colnames(dosage)))
print(paste0("number excluded (additive):",362 - ncol(dosage)))

covars <- covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4")]
covars["status"] = as.numeric(phe$V3)
covars <- covars[is.na(covars$status) == FALSE,]
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

if (dim(covars)[1] != dim(dosage)[1]) {
	dosage <- dosage[rownames(covars),]
}

results.gen <- vector("list", dim(dosage)[2])
names(results.gen) <- colnames(dosage)
results.add <- vector("list", dim(dosage)[2])
names(results.add) <- colnames(dosage)

if (run_dosage) {
#    for (gene in colnames(dosage[,1:30])) {
    for (gene in colnames(dosage)) {
        if (print_gene) {
            print(gene)
        }
    	covars["gcounts"] <- as.numeric(dosage[,gene])
    	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
    		       family="binomial", data=covars))
    	fit$deviance.resid <- NULL
    	results.add[[gene]] <- fit
    }
    
    saveRDS(results.add, paste0("output/reg_test/reg_results/", outname, "_add.rds"))
}

if (run_rounded) {

    if (test) {
        dosage <- readRDS("output/make_dosage_rds/test_100_rounded_remove.rds")
    } else {
        dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_rounded_remove.rds")
    }

    if (dim(covars)[1] != dim(dosage)[1]) {
    	dosage <- dosage[rownames(covars),]
    }

    # use same columns as before (based on number of nonzero entries)
    dosage <- dosage[,col_idx]
    print(paste0("number excluded (rounded additive):",362 - ncol(dosage)))
    #print(dosage)    
    #for (gene in colnames(dosage)) {
    #    print(length(unique(dosage[,gene])))
    #}
    if (rounded_add) {
#        for (gene in colnames(dosage[,1:30])) {
        for (gene in colnames(dosage)) {
            if (print_gene) {
                print(gene)
            }
            #print(paste0("gcounts rows: ", nrow(covars[,"gcounts"])))
            #print(paste0("dosage rows: ", nrow(as.numeric(dosage[,gene]))))
        	covars["gcounts"] <- as.numeric(dosage[,gene])
        
        	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
        		       family="binomial", data=covars))
        	fit$deviance.resid <- NULL
        	results.add[[gene]] <- fit

        }
        saveRDS(results.add, paste0("output/reg_test/reg_results/", outname, "_rounded_add.rds"))

    }
    if (rounded_factor) {
        # only include columns with more than 2 values (not including NA)
        unique_counts <- apply(dosage, 2, function(x)length(unique(x[!is.na(x)])))
        dosage <- dosage[,unique_counts > 2]
        print(paste0("number excluded (rounded factor):",362 - ncol(dosage)))

#        for (gene in colnames(dosage[,1:30])) {
        for (gene in colnames(dosage)) {
            if (print_gene) {
                print(gene) 
            }
        	covars["gcounts"] <- as.factor(dosage[,gene])
            n_levs <- nlevels(covars["gcounts"][which(!covars["gcounts"] == 0),])
           # print(covars["gcounts"][which(!covars["gcounts"] == 0),])
           # print(n_levs) 
            if (n_levs > 1) {
        	    fit <- summary(glm(status ~ as.factor(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
        	    	       family="binomial", data=covars))
        	    fit$deviance.resid <- NULL
        	    results.gen[[gene]] <- fit
            }
        }
        
        saveRDS(results.gen, paste0("output/reg_test/reg_results/", outname, "_rounded_gen.rds"))
    }
}

if (run_BMA) {

    # read in rounded dosage file
    if (test) {
        dosage <- readRDS("output/make_dosage_rds/test_100_rounded_remove.rds")
    } else {
        dosage <- readRDS("output/make_dosage_rds/ukb_hla_v2_rounded_remove.rds")
    }

    if (dim(covars)[1] != dim(dosage)[1]) {
    	dosage <- dosage[rownames(covars),]
    }

    #print(dosage)
#    BMA_haps_df <- as.data.frame(read.table("../notebooks/output/compare_cutoff_pvals/BMA_haps.csv", header=FALSE, sep = ",", row.names=1, stringsAsFactors=FALSE))

    # run on all haplotypes
    BMA_haps_df <- as.data.frame(read.table("../notebooks/output/compare_cutoff_pvals/BMA_haps_100.csv", header=FALSE, sep = ",", row.names=1, stringsAsFactors=FALSE))

    hap_names <-  as.character(BMA_haps_df[phe_name,])
    hap_names <- hap_names[!hap_names == ""]
    #print(BMA_haps_df)
    print(paste0("hap names: ", hap_names))

#    lines <- readLines("../notebooks/output/compare_cutoff_pvals/BMA_haps.csv")
#    for (i in 1:length(lines)) {
#        split_line <- strsplit(lines[i], ",")
#        temp_phe <- split_line[1]
#        print("temp_phe: ")
#        print(temp_phe)
#        print("alleles: ")
#        print(split_line[2:length(split_line)])
#    }

#    print(paste0("first lines: ",lines[0:5]))


    # choose the num_haps allelotypes with the smallest pvalues

    old_method <- FALSE

    if (old_method) {
        phe_results <- as.data.frame(read.table(paste0("output/process_output/", phe_name, "_processed_firth.txt"), header=FALSE, skip=1))

        phe_results <- phe_results[order(phe_results["V3"])[1:num_haps],]

   
        print(phe_results)
        print(phe_results["V1"])

        # get the names of the allelotypes we're using
        hap_names = phe_results[1:num_haps, "V1"]
        print(hap_names)
    }
    
    # subset dataframe to only the num_haps columns
    dosage <- dosage[, names(dosage) %in% hap_names]
    print(colnames(dosage))


    print("printed dosage")
    library("BMA")

    print("starting BMA")

    # run BMA analysis
    glm.out.FF <- bic.glm(dosage, covars$status, strict = FALSE, OR = 30,
        glm.family="binomial", factor.type=FALSE) 
    fit <- summary(glm.out.FF)
    print(paste0("postmean: ", glm.out.FF$postmean))
    print(paste0("postsd: ", glm.out.FF$postsd))

    # save output
    saveRDS(fit, paste0("output/reg_test/results_BMA/bma_", phe_name, "_",num_haps,"_round_all",".rds"))
    pdf(paste0("output/reg_test/plots_bma/bma_", phe_name, "_round_adjp_all",".pdf"))

    # plot image
    imageplot.bma(glm.out.FF)
    dev.off()
}

