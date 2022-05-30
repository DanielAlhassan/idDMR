#' Borrows DMRcate's 'KernelSums' code and applies a normalized Kernel-weight due to Alhassan et al (2022)
#' Uses an array-adaptive Kernel size
#' bw_opt() computes the median probe gaps per chromosome and uses this as the bandwidth in Kernel smoothing process
#' @param y refers to the observed chi-square values obtained from squaring limma's moderated t
#' @param x refers to the CpG site genomic position on the chromosome.
#' @param K refers to the Kernel which is Gaussian.
#' @noRd

bw_opt <- function(x) {
  bw <- median(na.omit(x)) # use median probe gaps.
  # Uneven distribution of probe spacing on each chromosome
  return(bw)
}

NormKernelSums <-
  function(x, y, K) {
    if (missing(K)) {
      K <- function(dx) {
        exp(-dx * dx / 2)
      } ## Gaussian
      attr(K, "support") <- 5
    }
    support <- attr(K, "support")

    stopifnot(all(diff(x) >= 0))

    # Bandwidth
    scale <-
      bw_opt(diff(x)[-1]) # bandwidth for each chromosome taken to be median probe gaps


    maxdiff <- support * scale
    triples <- SDeltas(x, maxdiff = maxdiff)
    # View(triples)

    i <- triples[, "i"]
    j <- triples[, "j"]
    dx <- triples[, "dx"]

    k <- K(dx / scale)

    # Normalize Kernel-weights
    df <- data.frame(i, j, dx, k)
    df2 <- df %>%
      dplyr::group_by(i) %>%
      dplyr::summarize(total = sum(k)) # total = Sum(k)

    df <- merge(x = df,
                y = df2,
                by = "i",
                all.x = TRUE)

    k <- df[, "k"] / df[, "total"] # k/sum(k) = w
    df$k <- k



    levels <- seq(length(x))
    ifac <- factor(i, levels = levels)
    jfac <- factor(j, levels = levels)

    itab <- table(ifac)
    jtab <- table(jfac)

    izero <- which(itab == 0)
    jzero <- which(jtab == 0)

    K0 <- K(0)

    SUM <-
      function(iterms, jterms) {
        S1 <- tapply(iterms, ifac, sum)
        S1[izero] <- 0
        S2 <- tapply(jterms, jfac, sum)
        S2[jzero] <- 0
        # browser()
        # tmp <- cbind(S1, S2); View(tmp)
        S1 + S2
      }

    # browser()
    swy <-
      K0 * y + SUM(k * y[j], k * y[i]) # normalized kernel-weighted statistic
    sw <- K0 + SUM(k, k)
    k2 <- k * k
    sww <- K0 * K0 + SUM(k2, k2)

    h_opt <- rep(scale, length(swy))

    out <- cbind(swy, sw, sww, h_opt)

    # print(out)
    return(out)
  }
