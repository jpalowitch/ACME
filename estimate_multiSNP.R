# Takes a SNP matrix, expression response, and covariates and 
# estimates the corresponding multivariate ACME model via MLE
#
# Arguments:
#  SNPs - a matrix of genotypes; columns are samples
#  expr - a vector of *ACME-transformed* expression values
#  cvrt - a matrix of covariates; columns are samples
#  method - the optim method to use; if NULL, default will be used.
#           Only "BFGS" currently implemented. See ?optim for details.
#
# Returns a list object with components:
#  beta0hat - estimated beta0
#  betahat - estimated SNP coefficients
#  gamhat - estimated cvrt coefficients
#  sigma2hat - estimated noise variance
#  df - degrees of freedom
#-------------------------------------------------------------------------------



estimate_multiSNP <- function (SNPs, expr, cvrt, method = "BFGS") {
  
  # Prepping data
  expr <- as.vector(expr)
  ncov <- nrow(cvrt); nsnp <- nrow(SNPs); n <- length(expr)
  X <- rbind(rep(1, n), SNPs)
  
  # RSS function
  loglik <- function(BetaGam) {
    Beta <- BetaGam[1:(nsnp + 1)]
    Gam <- BetaGam[(nsnp + 2):(nsnp + 1 + ncov)]
    Resid <- expr - log(crossprod(X, Beta)) - crossprod(cvrt, Gam)
    return(sum(Resid^2))
  }
  
  # Gradient function
  Dloglik <- function(BetaGam) {
    Beta <- BetaGam[1:(nsnp + 1)]
    Gam <- BetaGam[(nsnp + 2):(nsnp + 1 + ncov)]
    Resid <- expr - log(crossprod(X, Beta)) - crossprod(cvrt, Gam)
    DBeta <- -2 * crossprod(t(X) / as.vector(crossprod(X, Beta)), Resid)
    DGam <- -2 * cvrt %*% Resid
    return(c(DBeta, DGam))
  }
  
  # Finding start value for beta0
  beta0start <- exp(mean(expr))
  parStart <- c(beta0start, rep(0, nsnp + ncov))
  
  if (method == "BFGS") {
    estimates <- optim(parStart, loglik, gr = Dloglik, method = "BFGS")$par
  } else {
    estimates <- optim(parStart, loglik)$par
  }
  
  # Calculating RSS and sigmahat
  RSSfinal <- loglik(estimates)
  DF <- n - nsnp - ncov - 1
  
  return(list("beta0hat" = estimates[1],
              "betahat" = estimates[2:(nsnp + 1)],
              "gamhat" = estimates[(nsnp + 2):(nsnp + ncov)],
              "sigmahat" = RSSfinal / DF,
              "DF" = DF))
  
}


#-------------------------------------------------------------------------------
# Sample data and usage

if (FALSE) {
  
  n <- 1000; ncov <- 5; p <- 0.2
  beta0 <- 1000; beta1 <- 100; beta2 <- 200
  cvrt <- matrix(rnorm(ncov * n), nrow = ncov)
  SNPs <- matrix(sample(c(0:2), n * 2, replace = TRUE), nrow = 2)
  X <- rbind(rep(1, n), SNPs); Beta = c(beta0, beta1, beta2)
  trueGam <- rnorm(ncov)
  expr <- crossprod(X, Beta) * exp(crossprod(cvrt, trueGam) + rnorm(n)) - 1
  log_expr <- log(1 + expr)
  
  estimate_multiSNP(SNPs, log_expr, cvrt)
  
}
