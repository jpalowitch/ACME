# Make a toy eQTL data set for a single tissue
# n - number of samples (default = 100)
# ngene - number of genes (default = 1000)
# nsnp - number of snps (default = 10000)
# ncov - number of covariates (default = 10)
# p - reference allele probability (default = .8)
# saveDir - directory in which to save data files

make_eQTL_data <- function (n, ngene, nsnp, ncov, p, saveDir)