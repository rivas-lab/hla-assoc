args <- commandArgs(TRUE)
phe <- args[1]
outname <- args[2]

print("start")
options(warn=1)

test <- F
run_rounded <- F
run_dosage <- F
rounded_add <- F
rounded_factor <- F
run_BMA <- T
print_gene <- T
t0 <- proc.time()
#ukbb_files <- "$SCRATCH/ukbb_files/"
ukbb_files <- "/scratch/users/jolivier/ukbb_files/"
if (test) {
    dosage <- readRDS("test_100_remove.rds")
    covars <- readRDS("test_100_covar_16698_remove.rds")
#    outname <- "test_out_dosage"
} else {
    dosage <- readRDS("ukb_hla_v2_remove.rds")
    covars <- readRDS(paste0(ukbb_files,"ukb_hla_v2_covar_16698_remove.rds"))
#    outname <- "out_dosage"
}

phe <- as.data.frame(read.table(phe, header=FALSE, row.names=1))

# limit phe to individuals we're considering
phe <- phe[rownames(covars),]

print("done loading")

# use only columns with more than 5 nonzero entries
col_idx = colSums(dosage != 0) > 5
dosage <- dosage[,col_idx]
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
    
    saveRDS(results.add, paste0("reg_results/", outname, "_add.rds"))
}

if (run_rounded) {

    if (test) {
        dosage <- readRDS("test_100_rounded_remove.rds")
    } else {
        dosage <- readRDS(paste0(ukbb_files,"ukb_hla_v2_rounded_remove.rds"))
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
        saveRDS(results.add, paste0("reg_results/", outname, "_rounded_add.rds"))

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
        
        saveRDS(results.gen, paste0("reg_results/", outname, "_rounded_gen.rds"))
    }
}

if (run_BMA) {
    library("BMA")
    #print(dosage)
    print("starting BMA")
    dosage = cbind(dosage, covars)
#    dosage = subset(dosage, select = -c(status))
    #dosage = subset(dosage, select = c(A_101, B_801, C_701)) 
    dosage = dosage[,1:50]
    #print(dosage)

#    print(covars["status"])
#    y <- as.vector(covars["status"])
#    print(y)
    glm.out.FF <- bic.glm(dosage, covars$status, strict = FALSE, OR = 30,
        glm.family="binomial", factor.type=FALSE) 
    summary(glm.out.FF)
    print(glm.out.FF$probne0)
    #png("bicglm.png", type="quartz")
    imageplot.bma(glm.out.FF)
    #dev.off()
    write(ncol(dosage), file="timing.txt", append=TRUE)
}

#print(proc.time() - t0)
write(proc.time() - t0, file="timing.txt", append=TRUE)
#write(paste(ncol(dosage), proc.time, sep=","), file="timing.txt", append=TRUE)
