args = commandArgs(TRUE)
data <- args[1] # expression data file
pdata <- args[2] # phenotype data file
bioIDs <- args[3] # sample bioproject IDs file
s <- args[4] # column name of stress data in pdata
t <- args[5] # column name of tissue data
f <- args[6] # column name of family data
save <- args[7] # file save name 

# Load Necessary Packages
# BiocManager::install("sva") # install sva package
# BiocManager::install("limma") # install limma package
library(sva)
library(limma)
library(tidyverse)

# Load Datasets
# set working directory for easy access af data
#setwd(dir) 

# read in expression data
mdata <- read.csv(data, header=T) 
rownames(mdata) <- mdata[,1] ; mdata[,1] <- NULL # reset rownames to Orthogroup columns

# read in phenotype data
pheno <- read.csv(pdata, header=T) 
rownames(pheno) <- pheno[,1] # reset rownames to SRA accessions
 
# remove SRA accessions in pheno data that do not have expression data
pheno2 <- pheno[(rownames(pheno) %in% colnames(mdata)),] # SRAs in pheno that match mdata
pheno2$sample <- 1:nrow(pheno2) # create a sample column

# rename columns
pheno2 <- pheno2 %>% rename(stress=s, tissue=t, family=f)

# look at the number of samples in each factor
table(pheno2$tissue) ; table(pheno2$stress) ; table(pheno2$family)

# read in Bioproject ID information that corresponds to pheno2
bioproject <- read.delim(bioIDs, header=F, sep=",")
colnames(bioproject) <- c("sra", "bioproject") # rename columns

# combine bioproject and pheno2
pheno2 <- left_join(pheno2, bioproject, by="sra") # combine by matching sra columns
write.csv(pheno2, "pheno_data.csv")

# visualize expression and phenotype data
head(mdata[,1:5]) ; head(pheno2)

# Set Full and Null Models
# full model matrix - family, stress, tissue should be sig
mod <- model.matrix(~stress+tissue+family, data=pheno2)

# null model matrix (no adjustment variables are included)
null_mod <- model.matrix(~1, data=pheno2)

# expression data must be a matrix
mdata <- as.matrix(mdata)

# Perform SVASeq on Entire Dataset
# Estimate surrogate variables (SVs) using the two-step SVA method
svseq_obj <- svaseq(mdata, mod, null_mod, method = "two-step")
svseq_obj$sv[1:5, 1:5] # visualize surrogate variables

# plot of first 2 surrogate variables
png("SV_plot.png")
plot(svseq_obj$sv, pch=20, col="blue") ; dev.off()

# Adjust Data for SVs with Limma Package
## Remove effect of surrogate variables with limma package by fitting a linear model with surrogate variables included
# Full model with SVs
modSv <- cbind(mod, svseq_obj$sv)

# Null model with SVs as adjustment variables
mod0Sv <- cbind(null_mod, svseq_obj$sv) 

# linear model
fit <- lmFit(mdata, modSv) ; summary(fit)

# PCA
## Generate a clean matrix using a function by Andrew Jaffe
# This function removes the effects of SVs from our expression data
#'y' as the gene expresion matrix
#'mod' as the model matrix you sent to sva (the full model)
#'svs' as svobj$sv where svobj is the output from the sva function
cleanY = function(y, mod, svs) {
    X = cbind(mod, svs) # same as modSv
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

clean_data <- cleanY(mdata, mod, svseq_obj$sv)

# Function for plotting pca
plot_pca <- function(Legend, title){
  x %>% as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2, col=Legend)) + geom_point() + 
  ggtitle(title) +
  labs(x=paste("PC1: ",round(var_explained[1]*100,2),"%"),
       y=paste("PC2: ",round(var_explained[2]*100,2),"%"))
}
# Plot PCA of cleaned matrix
pca.res <- prcomp(t(clean_data), center = T, scale = T) # run PCA
write.csv(pca.res$x, "PCA_clean.csv") # save matrix of principle components (PCs)
var_explained <- pca.res$sdev^2/sum(pca.res$sdev^2) # explained variance for each PC
write.csv(var_explained, "PCA_var_clean.csv") # save matrix of explained variance per PC
x <- data.frame(pca.res$x) # create a dataframe of PCs
png("PCA_stress_clean.png")
plot_pca(pheno2$stress, "PCA on Clean Data (Stress)") ; dev.off()
png("PCA_tissue_clean.png")
plot_pca(pheno2$tissue, "PCA on Clean Data (Tissue)") ; dev.off()
png("PCA_family_clean.png")
plot_pca(pheno2$family, "PCA on Clean Data (Family)") ; dev.off()

# write clean data to file
write.csv(clean_data, paste(save, "_clean.csv", sep="")

# Adjusting for SVs using `f.pvalue` Method from SVA Package
# The f.pvalue function can be used to calculate parametric F-test p-values for each row of a data matrix
# The F-test compares the models mod and null_mod. They must be nested models, so all of the variables in null_mod must appear in mod.

# Calculate the F-test p-values for differential expression without adjusting for surrogate variables
pValues = f.pvalue(mdata,mod,null_mod)
qValues = p.adjust(pValues,method="BH")
write.csv(cbind(pValues, qValues), "Diff_expr_orthogroups_no_sv_adjustment.csv")

# Include the surrogate variables in both the null and full models to adjust for the surrogate variables by treating them as adjustment variables that must be included in both models. 
pValuesSv = f.pvalue(mdata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
write.csv(cbind(pValuesSv, qValuesSv), "Diff_expr_orthogroups_sv_adjusted.csv")

# Identify orthogroups (OGs) that are NOT diferentially expressed with respect to stress + tissue + family
length(pValuesSv[pValuesSv>0.05]) # 38 w/BH corrected p-value > 0.05
write.csv(pValuesSv[pValuesSv>0.05], "Not_diff_expr_orthogroups.csv") # corresponding OGs

# Correlation between SVs and Sample, Species, & BioProject
# Multiple linear regression for testing the association between SVs and 3 batch variables: bioproject, sample, and species
test <- lm(svseq_obj$sv ~ bioproject + sample + species, data = pheno2)
summary(test)

########################################################################
SVs <- svseq_obj$sv[,1:2] # create a matrix of SVs

# Simple Linear Regression to test the association between bioproject and SVs
test2 <- lm(SVs~bioproject, data = pheno2) ; summary(test2)

# Simple Linear Regression to test the association between sample and SVs
test3 <- lm(SVs~sample, data = pheno2) ; summary(test3)

# Simple Linear Regression to test the association between species and SVs
test4 <- lm(SVs~species, data = pheno2) ; summary(test4)

