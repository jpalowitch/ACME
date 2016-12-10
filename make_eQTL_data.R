# Make a toy eQTL data set for a single tissue
# n - number of samples
# ngene - number of genes (default = 1000)
# nsnp - number of snps (default = 10000)
# ncov - number of covariates (default = 10)
# p - reference allele probability (default = .2)
# saveDir - directory in which to save data files (default = cwd)

make_eQTL_data <- function (n, ngene = 1000, nsnp = 10000, ncov = 10, 
                            p = 0.2, saveDir = getwd()) { 
  
  if (ngene > nsnp)
    stop("must set ngene <= nsnp\n")
  
  # Making SNP data
  snps <- sample(c(0, 1, 2), n * nsnp, replace = TRUE,
                   prob = c((1 - p)^2, 2 * p * (1 - p), p^2))
  snps <- matrix(snps, nrow = n)
  
  # Making covariates and covariate effect
  cvrts <- matrix(rnorm(n * ncov), nrow = n)
  gam <- rnorm(ncov)
  cvrtEffect <- as.vector(cvrts %*% gam)
  
  # Making eQTL network
  causalSnps <- sample(1:nsnp, ngene, replace = FALSE)
  
  # Making effect sizes
  etas <- (exp(rnorm(ngene, 0, 1)) - 1) / 2
  
  # Mean function
  meanFun <- function (snpVec) {log(100) + log(1 + etas * snpVec)}
  
  # Making gene data
  genes <- snps[ , causalSnps]
  genes <- t(apply(genes, 1, meanFun)) + cvrtEffect +
    matrix(rnorm(n * ngene), nrow = n)
  genes <- exp(genes)
  
  return(list("snps" = snps, "genes" = genes, "cvrts" = cvrts))
  
}