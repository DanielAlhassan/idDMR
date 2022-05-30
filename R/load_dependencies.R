#' Load idDMR package dependencies
#' @description This function loads (and/or installs) dependencies needed
#' to use the idDMR package.
#'
#'
#' @usage load_dependencies()
#'
#' @export
#'
#' @author Daniel Alhassan
#'

load_dependencies = function() {

  packages = c("plyr", "BiocManager","dplyr")
  biocpackages = c("minfi", "limma", "ExperimentHub", "Gviz"
                   ,"edgeR","IlluminaHumanMethylation450kanno.ilmn12.hg19"
                   ,"IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  pkgs_to_install = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    } else{
      library(x, character.only = TRUE)
    }
  }
  biocpkgs_to_install = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
    } else{
      library(x, character.only = TRUE)
    }
  }

  ## Now load or install & load all
 lapply(packages, pkgs_to_install)
 lapply(biocpackages, biocpkgs_to_install)

 lapply(packages, require, character.only = TRUE)
 lapply(biocpackages, require, character.only = TRUE)

}


