#' Annotate Illumina CpGs with their chromosome position and test statistic based on \code{limma}.
#' Annotate a matrix or GenomicRatioSet of 450K or EPIC data with probe weights and chromosomal position.
#' @usage cpgsite.annotate(datatype = c("array", "sequencing"), object, what=c("Beta", "M"),
#' arraytype=c("EPIC", "450K"), analysis.type = "differential",
#'                                               design, contrasts = FALSE,
#'                                               cont.matrix = NULL, fdr = 0.05, coef, ...)
#' @param datatype Specify whether the data is of type "array" or "sequencing". Currently
#' only "array" works.
#' @param object A matrix-like or GRanges-like object with rows equal to the number of CpG sites interrogated and columns equal the number of samples
#' @param what Specify whether matrix-like or GRanges-like object contains \eqn{\beta} values or M-values
#' @param arraytype 450K or EPIC data. Only works for array \code{datatype}.
#' @param analysis.type Options include \code{differential} for \code{aaDMR()} to return DMRs;
#'                       \code{variability} to return VMRs. \code{variability} option yet to be implemented.
#' @param design Study design matrix. Identical context to differential analysis pipeline in \code{limma}. Must have an intercept if contrasts=FALSE. Applies only when analysis.type %in% c("differential").
#' @param fdr FDR cutoff (Benjamini and Hocheberg(1995)) for which CpG sites are individually declared significant. Used to index default thresholding in aaDMR(). \strong{Highly recommended as the primary thresholding parameter for calling DMRs}.
#'
#' @return A [CpGsiteAnnotated-class] class
#'
#' @export
#'


cpgsite.annotate <-
  function(datatype = c("array", "sequencing"),
           object,
           what = c("Beta", "M"),
           arraytype = c("EPIC", "450K"),
           analysis.type = "differential",
           design,
           contrasts = FALSE,
           cont.matrix = NULL,
           fdr = 0.05,
           coef,
           ...) {
    analysis.type <- match.arg(analysis.type)
    what <- match.arg(what)
    arraytype <- match.arg(arraytype)
    if (datatype == "array") {
      stopifnot(class(object)[1] %in% c("matrix", "GenomicRatioSet"))
      if (is(object, "matrix")) {
        if (arraytype == "450K") {
          grset <- makeGenomicRatioSetFromMatrix(
            mat = object,
            array = "IlluminaHumanMethylation450k",
            annotation = "ilmn12.hg19",
            what = what
          )
        }
        if (arraytype == "EPIC") {
          grset <- makeGenomicRatioSetFromMatrix(
            mat = object,
            array = "IlluminaHumanMethylationEPIC",
            annotation = "ilm10b4.hg19",
            what = what
          )
        }
      } else {
        grset <- object
      }
      object <- getM(grset)
      switch(analysis.type,
             differential = {
               stopifnot(is.matrix(design))
               if (!contrasts) {
                 stopifnot(colnames(design)[1] == "(Intercept)")
               } else {
                 stopifnot(!is.null(cont.matrix))
               }
               fit <- lmFit(object, design, ...)
               if (contrasts) {
                 stopifnot(coef %in% colnames(cont.matrix))
                 fit <- contrasts.fit(fit, cont.matrix)
               }
               fit <- eBayes(fit)
               tt <- topTable(fit, coef = coef, number = nrow(object))
               nsig <- sum(tt$adj.P.Val < fdr)
               if (nsig == 0) {
                 message(
                   "Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in aaDMR() to return DMRs, but be warned there is an increased risk of Type I errors."
                 )
               }
               if (nsig > 0 & nsig <= 100) {
                 message(
                   paste(
                     "Your contrast returned",
                     nsig,
                     "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."
                   )
                 )
               }
               if (nsig > 100) {
                 message(
                   paste(
                     "Your contrast returned",
                     nsig,
                     "individually significant probes. We recommend the default setting of pcutoff in aaDMR()."
                   )
                 )
               }
               betafit <- lmFit(ilogit2(object), design, ...)
               if (contrasts) {
                 betafit <- contrasts.fit(betafit, cont.matrix)
               }
               betafit <- eBayes(betafit)
               betatt <-
                 topTable(betafit, coef = coef, number = nrow(object))
               m <- match(rownames(tt), rownames(betatt))
               tt$diff <- betatt$logFC[m]
               m <- match(rownames(object), rownames(tt))
               tt <- tt[m,]
               anno <- getAnnotation(grset)
               stat <- tt$t
               annotated <-
                 GRanges(
                   as.character(anno$chr),
                   IRanges(anno$pos, anno$pos),
                   stat = stat,
                   diff = tt$diff,
                   ind.fdr = tt$adj.P.Val,
                   is.sig = tt$adj.P.Val < fdr
                 )
               names(annotated) <- rownames(tt)
             })
      annotated <- sort(annotated)
      return(new("CpGsiteAnnotated", ranges = annotated))
    }
    if (datatype == "sequencing") {
      stop("Sequencing mode is not functional at the moment.")
    } else {
      message("Error: datatype must be one of 'array' or 'sequencing'")
    }
  }
