# idDMR Package  <img src="images/logo.png" align="right" width="20%" height="20%" />
idDMR - Identify Differentially Methylated Regions for Microarray Data 🧬

<!-- badges: start -->
[![R](https://github.com/DanielAlhassan/idDMR/actions/workflows/r.yml/badge.svg)](https://github.com/DanielAlhassan/idDMR/actions/workflows/r.yml)
[![packageversion](https://img.shields.io/badge/Package%20version-0.5.0-orange.svg?style=flat-square)](commits/develop)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

<!-- badges: end -->


## Documentation


## Installation
<!--
You can install the release version from CRAN

``` r
install.packages("idDMR", dependencies=TRUE)
```

and the development version from GitHub
-->
``` r
install.packages("devtools")
devtools::install_github("DanielAlhassan/idDMR") 
```

## Usage
To use the **idDMR** package, start by loading dependencies which will install (where necessary) and load all dependent packages:

```r
library(idDMR)

load_dependencies()
```

To run site-level test, use the function **cpgsite.annotate()**
```r
myannotation <- cpgsite.annotate(datatype = "array",new_mval,what = "M",arraytype = "450K",
                                 analysis.type = "differential",design = design_mat,
                                 coef = 2, fdr = 0.05)
```

To 
```r
aadmr = aaDMR(myannotation,h = 1000,  min.cpgs = 2)
aadmr_df <- arrange(data.frame(extractRanges(aadmr, genome = "hg19")),seqnames)

```


## Packages
As well as the core idDMR, installing this package also installs a selection of other packages that you’re likely to use frequently. This includes packages for:


## Acknowledgments

<!--
## Citation
```r
@Package{,
  title = {idDMR: Identify Differentially Methylated Regions for Microarray Data},
  author = {Daniel Alhassan}, {Ebenezer Agbozo}
  year = {2022},
  note = {R package version 1.0.0},
  url = {https://github.com/DanielAlhassan/idDMR},
}

```
-->

## License
This package is free and open source software, licensed under GPL-3.
