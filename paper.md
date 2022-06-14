---
title: 'idDMR: Identify Differentially Methylated Regions for Microarray Data '
tags:
- R
- bioinformatics
- differentially methylated region
- kernel smoothing
- Illumina
date: "06 June 2022"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Daniel Alhassan
  corresponding: yes
  orcid: 0000-0003-0944-5818
  equal-contrib: yes
  affiliation: 1
- name: Gayla Olbricht
  orcid: 0000-0002-1213-2241
  equal-contrib: yes
  affiliation: 1
- name: Akim Adekpedjou
  orcid: 0000-0001-9584-4297
  equal-contrib: yes
  affiliation: 1
- name: Ebenezer Agbozo
  orcid: 0000-0002-2413-3815
  equal-contrib: yes
  affiliation: 2
bibliography: paper.bib
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
affiliations:
- name: Department of Mathematics and Statistics, Missouri University of Science and
    Technology, USA
  index: 1
- name: Department of Big Data Analytics and Methods of VideoAnalytics, Institute
    of               Radioelectronics and Information Technology, Ural Federal University
  index: 2
---

# Summary
Identifying and detecting differentially methylated regions(DMRs) in  the human genome has been an area of keen research for years now. Many methods and their associated `R` packages to solve problems related to methylated data. However, our knowledge on the methylation data particularly, the co-methylation structure between nearby CpG sites has increased. [@Sun2019b; @Sun2019] recently revealed that methylation patterns between neighboring sites is high for sites within a few hundred base pairs (bp) of each other. Even more recently in [@sun2022], it was revealed in a breast cancer study that the degree of co-methylation is different on different chromosomes. Owing to these recent findings, we developed the `idDMR`  R package to help the scientists in the field to implement the DMR detection method that uses an array-adaptive normalized kernel-weighted statistic developed by [Alhassan et al 2022].
`idDMR`[idDMR](https://github.com/DanielAlhassan/idDMR) is currently hosted in GitHub repository. Currently this package enhouses one DMR detection method in the function `aaDMR()`. This function uses the array-adaptive normalized kernel statistic in (\ref{A}). Some of its arguments, it accepts a `CpGsiteannotated` object, an adapted object in the `idDMR` package, from the originally created `cpgannotated` object by [@DMRcate2015]; an agglomerate parameter, g, set by the user that is used as the maximum bp distance in collapsing contiguous sites into regions. We describe the key functions in the package and how to use it in the Usage section below.



# Statement of need
Identifying differentially methylated regions is a common but challenging task the researcher faces owing to the complex nature and degree of co-methylation between neighboring CpG sites on a chromosome. Several efforts have been made, methods have been proposed and software packages written for these methods to make the users task easier. As our knowledge on the degree of co-methylation increases, better methods that utilize the new findings are proposed and hence software need to be created to make them available to the user to implement.

Some existing R packages in this area of research are the  `bumphunter` package [@bumphunter;@minfi], the `ChAMP` package with implements the `Probe Lasso` method to DMR identification
 
 As with the above R/Biconductor packages,  `idDMR` package was designed to allow the user to implement the methods discussed in [Alhassan et 2021]. At the heart of the detection procedure is the use of the normalized kernel-weighted statistic to incorporate co-methylation information at the CpG site-level. Our modeling framework borrows it's foundation work from DMRcate [@DMRcate2015].
 

### Mathematics
The estimator at the heart of our DMR detection procedure is the array-adaptive normalized kernel-weighted statistic $S(x_{i})$ below:

\begin{align}
S(x_{i}) = Y_{i} + \sum_{j \neq i}^{n}w_{j}(x_{i})Y_{j}
\label{A}
\end{align}


where $Y_{i}$ is the square of `limma`'s t-statistic [@limma] and $w_{j} = \dfrac{K\left( \dfrac{x_{j}-x_{i}}{h}\right)}{\displaystyle\sum_{j \neq i}^{n}K\left( \dfrac{x_{j}-x_{i}}{h}\right)}$.
See [Alhassan] for a more in dept discussion of the methods.

# Key Functions
The two core functions of the `idDMR` package are:
- `cpgsite.annotate` which runs CpG site-level test i.e. determine differentially methylated loci.

- `aadmr` which is responsible for identifying differentially methylated regions.


# Usage
To use the **idDMR** package, start by loading dependencies which will install (where necessary) and load all dependent packages:

```r
library(idDMR)

load_dependencies()
```

To run CpG site-level test (or determine differentially methylated loci), use the function `cpgsite.annotate()`
```r
myannotation <- cpgsite.annotate(datatype = "array", mval, what = "M", arraytype = "450K",
                                 analysis.type = "differential", design = design_mat,
                                 coef = 2, fdr = 0.05)
```

Next, use `aadmr()` function to identify differentially methylated regions.
```r
aadmr = aaDMR(myannotation, g = 1000,  min.cpgs = 2)

#extract the DMResults and output as a dataframe
aadmr_df <- arrange(data.frame(extractRanges(aadmr, genome = "hg19")), seqnames)

```

# Citations




# References
For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"
