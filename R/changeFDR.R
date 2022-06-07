#' Change the individual CpG FDR thresholding for a CpGsiteAnnotated object
#'
#' @description Takes a [CpGsiteAnnotated-class] object and a specified `0 < FDR< 1`, and re-indexes the object in order to call DMRs at the specified rate.
#'
#' @usage changeFDR(annot, FDR)
#'
#' @param annot A [CpGsiteAnnotated-class] object created from [cpgsite.annotate].
#' @param FDR The desired individual CpG FDR, which will index the rate at which DMRs are called.
#' @export
#'
#' @details The number of CpG sites called as significant by this function will set the post-smoothing threshold for DMR constituents in [aaDMR]
#' @return A re-indexed [CpGsiteAnnotated-class] object.


changeFDR <- function (annot, FDR)
{
  if(!is(annot, "CpGannotated")){
    stop("Error: annot is not a CpGsiteAnnotated object. Please create one with cpgsite.annotate()")
  }
  if(FDR <=0 | FDR >=1){
    stop("Error: please enter an appropriate FDR value, 0 < FDR < 1.")
  }
  annot@ranges$is.sig <- annot@ranges$ind.fdr < FDR
  cat(paste0("Threshold is now set at FDR=", FDR, ", resulting in ",
             sum(annot@ranges$is.sig), " significantly differential CpGs."))
  annot

}
