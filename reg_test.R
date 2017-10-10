args <- commandArgs(TRUE)
phe <- args[1]

test <- F
run_rounded <- T
run_dosage <- F
t0 <- proc.time()

if (test) {
    dosage <- readRDS("test_100_remove.rds")
    covars <- readRDS("test_100_covar_16698_remove.rds")
    outname <- "test_out_dosage"
} else {
    dosage <- readRDS("ukb_hla_v2_remove.rds")
    covars <- readRDS("ukb_hla_v2_covar_16698_remove.rds")                          
    outname <- "out_dosage"
}

phe <- as.data.frame(read.table(phe, header=FALSE, row.names=1))

print(nrow(dosage))
print(nrow(covars))
phe <- phe[rownames(covars),]
print(nrow(phe))

print("done loading")

# use only columns with more than 5 nonzero entries
#print(colSums(dosage != 0 & dosage != NA) > 5)
dosage <- dosage[,colSums(dosage != 0) > 5]

print(paste0("rows in dosage: ", nrow(dosage)))
print(paste0("rows in covars: ", nrow(covars)))


print("number excluded:")
print(362 - ncol(dosage))
#print(colSums(dosage != 0))

#rounded_dosage <- readRDS("ukb_hla_v2_rounded_remove.rds")
#for (gene in colnames(dosage)) {
#    print(unique(dosage[,gene]))
###    print(nnzero(dosage[,gene]))
##    print(colSums(dosage != 0))
#}

#print(phe[1:10,])
#print(dosage[1:10,])
#print(rounded_dosage[1:10,])
#print(covars[1:10,])

covars <- covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4")]
covars["status"] = as.numeric(phe$V3)
covars <- covars[is.na(covars$status) == FALSE,]
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

print(paste0("dim covars: ", nrow(covars)))

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

            print(gene)
    	covars["gcounts"] <- as.numeric(dosage[,gene])
    	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
    		       family="binomial", data=covars))
    	fit$deviance.resid <- NULL
    	results.add[[gene]] <- fit
    }
    
    #print(results.add)
    saveRDS(results.add, paste(outname, "_add.rds", sep=""))
    #save(results.add, paste(outname, "_add.txt", sep=""))
}

if (run_rounded) {

    if (test) {
        dosage <- readRDS("test_100_rounded_remove.rds")
    } else {
        dosage <- readRDS("ukb_hla_v2_rounded_remove.rds")
    }

    
    if (dim(covars)[1] != dim(dosage)[1]) {
    	dosage <- dosage[rownames(covars),]
    }


    
    #for (gene in colnames(dosage)) {
    #    print(length(unique(dosage[,gene])))
    #}

    unique_counts <- apply(dosage, 2, function(x)length(unique(x[!is.na(x)])))
    
    dosage <- dosage[,unique_counts > 2]
    print(paste0("rows in dosage: ", nrow(dosage)))
    print("number excluded:")
    print(362 - ncol(dosage))

    for (gene in colnames(dosage)) {
        print(gene)
        #print(paste0("gcounts rows: ", nrow(covars[,"gcounts"])))
        #print(paste0("dosage rows: ", nrow(as.numeric(dosage[,gene]))))
        print(length(as.numeric(dosage[,gene])))
        print(nrow(covars))
    	covars["gcounts"] <- as.numeric(dosage[,gene])
        print("problem here?")
    
    	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
    		       family="binomial", data=covars))
    	fit$deviance.resid <- NULL
    	results.add[[gene]] <- fit
    
    	covars["gcounts"] <- as.factor(dosage[,gene])
    #        print(covars["gcounts"])
            #print(covars["gcounts"][which(!covars["gcounts"] == 0),])
            n_levs <- nlevels(covars["gcounts"][which(!covars["gcounts"] == 0),])
            #print("n_levs:")
            #print(n_levs)
            #print(n_levs)
    
            if (n_levs > 1) {
    	    fit <- summary(glm(status ~ as.factor(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
    	    	       family="binomial", data=covars))
    	    fit$deviance.resid <- NULL
    	    results.gen[[gene]] <- fit
            }
    }
    
    saveRDS(results.gen, paste(outname, "_rounded_gen.rds", sep=""))
    saveRDS(results.add, paste(outname, "_rounded_add.rds", sep=""))
}

print(proc.time() - t0)
