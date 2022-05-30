#' S4 objects summarizing individual CpG sites after fitting to model.
#' @description An S4 class that stores output from [cpgsite.annotate()]
#' @slot ranges A GRanges object, containing CpG-level information to be passed to \link{aaDMR}. Mcols of this object include:
#'      \itemize{
#'        \item{stat:} {Per-CpG test statistic; t from limma.}
#'        \item{diff:} {Methylation difference/coefficient. In beta space for [cpgsite.annotate] output }
#'        \item{ind.fdr:} {False discovery rate calculated on individual CpG sites.}
#'        \item{is.sig:} {Logical determining whether a CpG site is individually significant or not.}
#' }
#'
#' @method \code{cpgsite.annotate} \code{show}




setClass(
  "CpGsiteAnnotated",
  slots = c(ranges = "GRanges"),
  prototype = list(ranges = GenomicRanges::GRanges())
)


#' Storage object for DMResults
#' @description An S4 class that stores DMR information as output from \link{aaDMR}.
#' @slot DMResults This class has eight slots, summarizing DMR information to be passed to \code{\link{extractRanges}}:
#' \itemize{
#'      \item{coord:} {DMR coordinates in UCSC style.}
#'      \item{no.cpgs:} {Number of consecutive CpG sites of DMR.}
#'      \item{min_smoothed_fdr:} {Minimum FDR of the smoothed statistic.}
#'      \item{Stouffer:} {p-value based on Stouffer's summary transform of the \strong{individual} CpG FDRs.}
#'      \item{Fisher:} {p-value based on Fisher combined probability transform of the \strong{individual} CpG FDRs}
#'      \item{maxdiff:} {Maximum differential within the DMR.}
#'      \item{meandiff:} {Mean differential across the DMR.}
#'      }
#' @method \code{DMResults}  \code{show}
#'
#'
setClass(
  "DMResults",
  slots = c(
    coord = "character",
    no.cpgs = "integer",
    min_smoothed_fdr = "numeric",
    Stouffer = "numeric",
    Fisher = "numeric",
    maxdiff = "numeric",
    meandiff = "numeric"
  ),
  prototype = list(
    coord = character(),
    no.cpgs = integer(),
    min_smoothed_fdr = numeric(),
    Stouffer = numeric(),
    Fisher = numeric(),
    maxdiff = numeric(),
    meandiff = numeric()
  )
)
