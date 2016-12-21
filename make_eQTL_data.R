# Make a toy eQTL data set for a single tissue
# dependencies: filematrix
# 
# Arguments:
#   n - number of samples
#   ngene - number of genes (default = 1000)
#   nsnp - number of snps (default = 10000)
#   ncov - number of covariates (default = 10)
#   p - reference allele probability (default = .2)
#   saveDir - directory in which to save data files (default = cwd)
#             (function will make a subdirectory named "filematrices")
#   returnData - whether or not to return data in a list (default = FALSE)
#                if TRUE, data is returned in a list with components:
#                  'snps'    - toy SNP data
#                  'genes'   - toy gene data
#                  'cvrts'   - toy covariate data
#                  'snploc'  - dummy SNP locations
#                  'geneloc' - dummy gene locations
#   fmat - if TRUE, will save data to filematrices. Otherwise, saves to 
#          text files. (default = FALSE)
#   verbose - whether or not to print SNP filematrix writing progress
#             (default = FALSE)
#-------------------------------------------------------------------------------

require(filematrix)

make_eQTL_data <- function (n, ngene = 1000, nsnp = 10000, ncov = 10, p = 0.2, 
                            saveDir = getwd(), returnData = FALSE, fmat = FALSE,
                            verbose = FALSE) { 
  
  if (ngene > nsnp)
    stop("must set ngene <= nsnp\n")
  
  # Making SNP data
  snps <- sample(c(0, 1, 2), n * nsnp, replace = TRUE,
                   prob = c((1 - p)^2, 2 * p * (1 - p), p^2))
  snps <- matrix(snps, ncol = n)
  
  # Making covariates and covariate effect
  cvrts <- matrix(rnorm(n * ncov), ncol = n)
  gam <- rnorm(ncov)
  cvrtEffect <- as.vector(crossprod(cvrts, gam))
  
  # Making eQTL network
  causalSnps <- sample(1:nsnp, ngene, replace = FALSE)
  
  # Making effect sizes
  etas <- (exp(rnorm(ngene, 0, 1)) - 1) / 2
  
  # Mean function
  meanFun <- function (snpVec) {log(100) + log(1 + etas * snpVec)}
  
  # Making gene data
  genes <- snps[causalSnps, ]
  genes <- apply(genes, 2, meanFun) + cvrtEffect +
    matrix(rnorm(n * ngene), ncol = n)
  genes <- exp(genes)
  
  # Naming gene/snp rows and making gene/snp loc dummy dataframes
  rownames(genes) <- 1:ngene
  rownames(snps) <- 1:nsnp
  rownames(cvrts) <- 1:ncov
  colnames(genes) <- colnames(snps) <- colnames(cvrts) <- 1:n
  geneloc <- data.frame("geneid" = rownames(genes),
                        "chrm_probe" = rep(1, ngene),
                        "s1" = c(1:ngene) * 10, "s2" = c(1:ngene) * 10)
  snpsloc <- data.frame("SNP" = rownames(snps),
                        "chrm_snp" = rep(1, nsnp),
                        "pos" = 1:nsnp)
  
  # Filtering gene data
  rpkm_counts <- rowSums(genes > 0.1)
  usable_genes <- rownames(genes)[rpkm_counts >= 10]
  gene_t1 <- genes[usable_genes, ]
  log_reads <- log(1 + genes)
  
  if (!fmat) {
    
    # Making directory for txtfiles and going there
    if (!dir.exists(file.path(saveDir, "txtfiles"))) {
      dir.create(file.path(saveDir, "txtfiles"))
    }
    
    write.table(geneloc, file = file.path(saveDir, "txtfiles", "gene_loc.txt"),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
    write.table(snpsloc, file = file.path(saveDir, "txtfiles", "snps_loc.txt"),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
    write.table(genes, file = file.path(saveDir, "txtfiles", "gene.txt"),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
    write.table(snps, file = file.path(saveDir, "txtfiles", "snps.txt"),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
    write.table(cvrts, file = file.path(saveDir, "txtfiles", "cvrt.txt"),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
  } else {
    
    # get single number locations: gene_loc, snps_loc
    {
      chrset = unique(c(unique(geneloc$chrm_probe), unique(snpsloc$chrm_probe)));
      geneloc$commonchr = match(geneloc$chrm_probe,   chrset, nomatch = -1L);
      snpsloc$commonchr = match(snpsloc$chrm_snp, chrset, nomatch = -2L);
      
      maxpos = 1e9;
      # maxpos = max( diff(range(geneloc$s1)), diff(range(snpsloc$pos)) ) + 1;
      
      geneloc$singlepos = geneloc$commonchr*maxpos + geneloc$s1;
      snpsloc$singlepos = snpsloc$commonchr*maxpos + snpsloc$pos;
      
      mch = match(rownames(log_reads), geneloc$geneid, nomatch = 0L)
      stopifnot(all(mch != 0L))
      gene_loc = geneloc$singlepos[mch];
      rm(mch);
      
      mch = match(rownames(snps), snpsloc$SNP, nomatch = 0L)
      stopifnot(all(mch != 0L))
      snps_loc = snpsloc$singlepos[mch];
      rm(mch);
      
      rm(maxpos, chrset);
    }
    
    # reorder genes / SNPs
    if( is.unsorted(gene_loc) ) {
      ord = sort.list( gene_loc );
      log_reads = log_reads[ord, ];
      gene_loc = gene_loc[ord];
      rm(ord);
    }
    
    if( is.unsorted(snps_loc) ) {
      ord = sort.list( snps_loc );
      snps = snps[ord, ];
      snps_loc = snps_loc[ord];
      rm(ord);
    }
    
    # Making directory for filematrices and going there
    if (!dir.exists(file.path(saveDir, "filematrices"))) {
      dir.create(file.path(saveDir, "filematrices"))
    }
    setwd(file.path(saveDir, "filematrices"))
    
    # Save gene_loc, snps_loc
    fm = fm.create.from.matrix('gene_loc', gene_loc)
    close(fm);
    fm = fm.create.from.matrix('snps_loc', snps_loc)
    close(fm);
    
    
    # Save SNPs in a filematrix
    fm = fm.create('snps', nrow = ncol(snps), ncol = nrow(snps), type = 'integer', size = 2);
    
    step1 = 100000;
    mm = nrow(snps);
    nsteps = ceiling(mm/step1);
    for( part in 1:nsteps ) { # part = 1
      if (verbose)
        cat( part, 'of', nsteps, '\n');
      fr = (part-1)*step1 + 1;
      to = min(part*step1, mm);
      
      fm[,fr:to] = t(snps[fr:to,])*1000;
      gc();
    }
    rm(part, step1, mm, nsteps, fr, to);
    
    colnames(fm) = rownames(snps);
    rownames(fm) = colnames(snps);
    close(fm)
    
    # Save genes 
    fm = fm.create.from.matrix('gene', t(log_reads))
    close(fm)
    
    # Save cvrt 
    fm = fm.create.from.matrix('cvrt', t(cvrts))
    close(fm)
    
    # Re-setting working directory
    setwd(saveDir)
  }
  
  if (returnData) {
    return(list("snps" = snps, "gene" = genes, "cvrt" = cvrts,
                "gene_loc" = genelocs, "snps_loc" = snplocs))
  }
  
}