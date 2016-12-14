# Takes a SNP matrix, expression response, and covariates and 
# estimates the corresponding multivariate ACME model via MLE






estimate_multiSNP <- function (SNPs, expr, cvrt) {
  
  expr <- as.vector(expr)
  
  
  
  
}


# Sample data and usage
n <- 100; ncov <- 5; p <- 0.2
beta0 <- 10000; beta1 <- 100; beta2 <- 200
cvrt <- matrix(rnorm(ncov * n), nrow = ncov)
SNPs <- matrix(sample(c(0:2), n * 2, replace = TRUE), nrow = 2)
X <- rbind(rep(1, n), SNPs); Beta = c(beta0, beta1, beta2)
log_expr <- log(crossprod(X, Beta)) + crossprod(cvrt, rnorm(ncov)) + rnorm(n)
