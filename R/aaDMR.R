#' DMR detection
#' @description The main function of this package. Computes a normalized kernel-weight statistic,
#'  \code{s(x)}, against a null comparison to identify significantly differentially (or variable) methylated regions.
#'
#' \deqn{ S(x_{i}) = Y_{i} + \sum_{\displaystyle\substack{j = 1\\ j \neq i}}^{n}w_{j}(x_{i})Y_{j}
#' } where
#' \eqn{ w_{j} = \dfrac{K\left( \frac{x_{j}-x_{i}}{h}\right)}{\sum_{\displaystyle\substack{j =1\\ j \neq i}}^{n}K\left( \frac{x_{j}-x_{i}}{h}\right)}}
#'
#' @usage aaDMR(object, h = 1000, pcutoff = "fdr", betacutoff = NULL, min.cpgs = 2)
#'
#' @param object A \link{CpGsiteAnnotated} object created from \link{cpgsite.annotate}.
#' @param h An agglomerate parameter. Any sig. cpgs within \code{h} distance are collapsed to form a DMR.
#' @param pcutoff A threshold to determine DMRs. Default is highly recommended.
#' @param min.cpgs The minimum number of consecutive CpGs constituting a DMR.
#'
#' @export
#'
#' @details The value of h should be chosen with care. We recommend using a maximum value of 1000bp
#' as the agglomerate parameter. If  `h` is too small DMRs might have refers CpGs and hence DMRs will
#' span shorter gene loci. Also if `h` is too large, DMRs will span larger gene loci with a increased risk of Type I errors.

aaDMR <- function(object,
                  h = 1000,
                  # h
                  pcutoff = "fdr",
                  betacutoff = NULL,
                  min.cpgs = 2) {
  ## Arguments
  stopifnot(is(object, "CpGsiteAnnotated"))
  stopifnot(h >= 1)
  stopifnot(pcutoff == "fdr" | (0 <= pcutoff & pcutoff <= 1))


  ## Modified 'object'
  object <- data.frame(
    ID = names(object@ranges),
    weights = abs(object@ranges$stat),
    CHR = seqnames(object@ranges),
    pos = start(object@ranges),
    diff = object@ranges$diff,
    indfdr = object@ranges$ind.fdr,
    is.sig = object@ranges$is.sig
  )


  # Order by position
  object <- object[order(object$CHR, object$pos),]


  ## Kernel (chi-squared) test via 'ParallelFits'
  lag <- h
  chr.unique <- unique(c(as.character(object$CHR)))
  fitted <-
    lapply(chr.unique,
           ParallelFits,
           object = object,
           h = h)
  object <- rbind.fill(fitted)

  ## FDR stuff
  object$fdr <- p.adjust(object$raw, method = "BH")
  if (pcutoff == "fdr") {
    nsig <- sum(object$is.sig)
    if (nsig == 0) {
      txt <-
        "The FDR you specified in cpgsite.annotate() returned no significant CpGs, hence there are no DMRs.\n    Try specifying a value of 'pcutoff' in aaDMR() and/or increasing 'fdr' in cpgsite.annotate()."
      stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
    }
    pcutoff <- sort(object$fdr)[nsig]
  }
  object$sig <- (object$fdr <= pcutoff)
  if (nrow(object) == 0) {
    txt <-
      "No signficant regions found. Try increasing the value of\n    'pcutoff' in aaDMR() and/or 'fdr' in cpg.annotate()."
    stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
  }

  ## Segmentation
  message("Identifying regions ...")
  # Define jump.k
  # K = number of significant CpGs
  # k = K - 1
  chr.N <- as.character(object$CHR)
  pos.N <- object$pos
  sig.N <- object$sig
  N <- length(sig.N)
  n.K <- which(sig.N)
  K <- length(n.K)
  stopifnot(K >= 2)
  pos.K <- pos.N[n.K]
  chr.K <- chr.N[n.K]
  jump_chr.k <- (chr.K[-1] != chr.K[-K])
  jump_pos.k <- (diff(pos.K) > lag)
  jump.k <- (jump_chr.k | jump_pos.k)
  # Segment using jump.k
  ksegments.A2 <- Segment(jump.k)
  A <- nrow(ksegments.A2)
  # Extract start/end indices
  kstart.A <- ksegments.A2[, "start"]
  kend.A <- ksegments.A2[, "end"]

  ## Regionwise stats
  # Positions
  realpos.K <- pos.K


  # Per-DMR: Coordinates
  start.A <- realpos.K[kstart.A]
  end.A <- realpos.K[kend.A]
  chr.A <- chr.K[kstart.A]
  stopifnot(all(chr.K[kend.A] == chr.A))
  fmt <- "%s:%1d-%1d"
  coord.A <- sprintf(fmt, chr.A, start.A, end.A)

  # Region factor
  nstart.A <- n.K[kstart.A]
  nend.A <- n.K[kend.A]
  width.A <- nend.A + 1 - nstart.A
  a.Z <- rep(seq(A), width.A) # a.Z
  fn <-
    function(a) {
      seq(from = nstart.A[a], to = nend.A[a])
    }
  l.listA <- lapply(seq(A), fn)
  n.Z <- unlist(l.listA)
  region.N <- rep(NA_integer_, N)
  region.N[n.Z] <- a.Z
  levels <- seq(A)
  region.N <- factor(region.N, levels = levels)

  # Per-DMR: Number of CpGs
  no_cpg.A <- c(table(region.N))

  # Function to do regionwise summaries
  REGIONSTAT <-
    function(field,
             fn) {
      x.N <- object[[field]]
      x.R <- tapply(x.N, region.N, fn)
      c(x.R)
    }
  # results <- region-wise stats
  fn_Stouffer <- function(x)
    pnorm(sum(qnorm(x)) / sqrt(length(x)))
  fn_Fisher <-
    function(x)
      pchisq((sum(log(x)) * -2), df = length(x) * 2, lower.tail = FALSE)
  fn_max <- function(x)
    x[which.max(abs(x))]
  results <-
    data.frame(
      coord = coord.A,
      no.cpgs = no_cpg.A,
      min_smoothed_fdr = REGIONSTAT("fdr", min),
      Stouffer = REGIONSTAT("indfdr", fn_Stouffer),
      Fisher = REGIONSTAT("indfdr", fn_Fisher),
      maxdiff = REGIONSTAT("diff", fn_max),
      meandiff = REGIONSTAT("diff", mean),
      row.names = seq(A),
      stringsAsFactors = FALSE
    )

  # Order and filter DMRs

  keep <- (results$no.cpgs >= min.cpgs)
  results <- results[keep,]
  if (!(is.null(betacutoff))) {
    message("Warning: betacutoff only meaningful for Illumina array data.")
    keep <- (abs(results$meandiff) > betacutoff)
    results <- results[keep,]
  }
  o <- order(results$Fisher,-results$no.cpgs)
  results <- results[o,]
  message("Fini!")
  return(
    new(
      "DMResults",
      coord = results$coord,
      no.cpgs = results$no.cpgs,
      min_smoothed_fdr = results$min_smoothed_fdr,
      Stouffer = results$Stouffer,
      Fisher = results$Fisher,
      maxdiff = results$maxdiff,
      meandiff = results$meandiff
    )
  )
}
