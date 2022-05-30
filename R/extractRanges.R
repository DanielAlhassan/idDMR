#' Create a GRanges object from \code{\link{aaDMR}} output.
#' @description Takes a DMResults object and produces the corresponding GRanges object.
#' @usage extractRanges(dmroutput, genome = c("hg19", "hg38", "mm10"))
#' @author Daniel Alhassan
#' @export
extractRanges <-
  function(dmroutput,
           genome = c("hg19", "hg38", "mm10")) {
    genome <- match.arg(genome)
    if (!is(dmroutput, "DMResults")) {
      stop("Error: dmroutput is not a DMResults object. Please create one with aaDMR().")
    }
    coords <- extractCoords(dmroutput@coord)
    coords <- cbind(
      coords,
      dmroutput@no.cpgs,
      dmroutput@min_smoothed_fdr,
      dmroutput@Stouffer,
      dmroutput@Fisher,
      dmroutput@maxdiff,
      dmroutput@meandiff
    )
    coords$chromStart <- as.integer(as.character(coords$chromStart))
    coords$chromEnd <- as.integer(as.character(coords$chromEnd))
    ranges <-
      makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
    eh <- ExperimentHub()
    switch(genome,
      hg19 = {
        grt <- eh[["EH3132"]]
      },
      hg38 = {
        grt <- eh[["EH3134"]]
      },
      mm10 = {
        grt <- eh[["EH3136"]]
      }
    )
    genesidx <- as.data.frame(findOverlaps(ranges, grt))
    genesover <-
      tapply(genesidx$subjectHits, genesidx$queryHits, function(x) {
        grt$symbol[x]
      })
    op.A <- sapply(genesover, function(l) {
      paste(l, collapse = ", ")
    })
    name.A <- names(genesover)
    m.A <- as.numeric(name.A)
    M <- length(ranges)
    overlapping.genes <- rep(NA_character_, M)
    overlapping.genes[m.A] <- op.A
    ranges$overlapping.genes <- overlapping.genes
    colnames(values(ranges)) <-
      sub("dmroutput@", "", colnames(values(ranges)))
    ranges
  }
