# Description: Obtain a list of SRA samples with expression data
args = commandArgs(TRUE)
data <- args[1] # expression data file
pdata <- args[2] # phenotype data file

# read in expression data
mdata <- read.csv(data, header=T) 
rownames(mdata) <- mdata[,1] ; mdata[,1] <- NULL # reset rownames to Orthogroup columns

# read in phenotype data
pheno <- read.csv(pdata, header=T) 
 
# remove SRA accessions in pheno data that do not have expression data
pheno2 <- pheno[(pheno$sra %in% colnames(mdata)),] # SRAs in pheno that match mdata

# save list of sras with expression data
write.table(pheno2$sra, "sra.txt", sep=",", row.names=F, col.names=F)