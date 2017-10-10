args <- commandArgs(TRUE)
phe <- args[1]

test <- TRUE

if (test) {
    dosage <- readRDS("test_100_remove.rds")
    rounded_dosage <- readRDS("test_100_rounded_remove.rds")
    covars <- readRDS("test_100_covar_16698_remove.rds")
    outname <- "test_out_dosage"
} else {
    dosage <- readRDS("ukb_hla_v2_remove.rds")
    rounded_dosage <- readRDS("ukb_hla_v2_rounded_remove.rds")
    covars <- readRDS("ukb_hla_v2_covar_16698_remove.rds")                          outname <- "out_dosage"
}

phe <- as.data.frame(read.table(phe, header=FALSE, row.names=1))

print(nrow(dosage))
print(nrow(rounded_dosage))
print(nrow(covars))
phe <- phe[rownames(covars),]
print(nrow(phe))

print("done loading")

#print(phe[1:10,])
#print(dosage[1:10,])
#print(rounded_dosage[1:10,])
#print(covars[1:10,])

covars <- covars[, c("age", "sex", "Array", "PC1", "PC2", "PC3", "PC4")]
covars["status"] = as.numeric(phe$V3)
covars <- covars[is.na(covars$status) == FALSE,]
covars[covars$status == 1, "status"] = rep(0, sum(covars$status == 1))
covars[covars$status == 2, "status"] = rep(1, sum(covars$status == 2))

if (dim(covars)[1] != dim(dosage)[1]) {
	dosage <- dosage[rownames(covars),]
        rounded_dosage <- dosage[rownames(covars),]
}

results.gen <- vector("list", dim(dosage)[2])
names(results.gen) <- colnames(dosage)
results.add <- vector("list", dim(dosage)[2])
names(results.add) <- colnames(dosage)

for (gene in colnames(dosage)) {
	covars["gcounts"] <- as.numeric(dosage[,gene])
	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
		       family="binomial", data=covars))
	fit$deviance.resid <- NULL
	results.add[[gene]] <- fit
}

#print(results.add)
saveRDS(results.add, paste(outname, "_add.rds", sep=""))
#save(results.add, paste(outname, "_add.txt", sep=""))


for (gene in colnames(rounded_dosage)) {
	covars["gcounts"] <- as.numeric(rounded_dosage[,gene])

	fit <- summary(glm(status ~ as.numeric(gcounts) + age + sex + Array + PC1 + PC2 + PC3 + PC4,
		       family="binomial", data=covars))
	fit$deviance.resid <- NULL
	results.add[[gene]] <- fit

	covars["gcounts"] <- as.factor(rounded_dosage[,gene])
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
