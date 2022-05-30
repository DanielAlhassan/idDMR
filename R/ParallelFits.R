#' Adapts DMRcate's fitParallel() function
#' Uses the array-adaptive normalized kernel-weight proposed by Alhassan et al (2022)
#' @param chr refers to the chromosome-type eg. chr1, chr2 etc.
#' @param h refers to agglomerate parameter. Any sig. cpgs with 'h' distance are collapsed to form a DMR
#' @noRd


ParallelFits <- function(chr, object, h) {
  message(paste("Array-adaptive Fitting for ", chr, "...", sep = ""))
  chrIndex <- object$CHR %in% chr
  chromosome <- object[chrIndex,]

  pos <- chromosome$pos


  lag <- h
  beta <- chromosome$weights
  df <- 1
  X2 <- beta^2
  j <- NormKernelTest(pos = pos, X2 = X2, df = df)
  pvalue <- j[, "pval"]

  chromosome$raw <- pvalue
  chromosome$h_opt <- j[, "h_opt"]

  # print(chromosome)
  return(chromosome)
}
