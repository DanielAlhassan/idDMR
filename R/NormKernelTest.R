#' Adapts DMRcate's 'KernelTest' function
#' Uses array-adaptive Kernel size in testing
#' @param pos refers to the CpG site genomic position on the chromosome
#' @param X2 refers to the observed chi-square values obtained from squaring limma's moderated t
#' @param df refers to the degree of freedom for X2 which is 1.
#' @keywords Kernel testing, Cochrane-Satterthwaite approximation
#' @noRd


NormKernelTest <-
  function(pos, # x values
           X2, # Chi-square values, each with 'df' degrees of freedom
           df = 1) {
    # Local kernel sums
    j <- NormKernelSums(x = pos, y = X2)
    swy <- j[, "swy"]
    sw <- j[, "sw"]
    sww <- j[, "sww"]
    h_opt <- j[, "h_opt"]

    # Cochrane-Satterthwaite approximation
    p <- df * sww / sw # typically between 0 and 1
    q <- sw * sw / sww # df for

    outX2 <-
      swy / p # Approx. chi-square with q degrees of freedom


    pval <- pchisq(outX2, df = q, lower.tail = FALSE)

    return(cbind(pval, h_opt))
  }
