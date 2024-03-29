---
title: 'idDMR: Identify Differentially Methylated Regions for Microarray Data '
tags:
- R
- bioinformatics
- differentially methylated region
- kernel smoothing
- Illumina
date: "September 17, 2022"
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
Identifying and detecting differentially methylated regions (DMRs) in  the human genome has been an area of keen research interest for years now. Many methods and their associated `R` packages have been proposed to solve problems related to DNA methylation data.  A few common existing R packages in this area of research are the  `bumphunter` package [@bumphunter;@minfi], the `ChAMP` package with implements the `Probe Lasso` method to DMR identification and the `DMRcate` package [@DMRcate2015]. \par
Our knowledge on the co-methylation structure between CpG sites has increased. [@Sun2019] recently revealed that co-methylation patterns within normal tissues  were as short as a few hundred base pairs (bp) of each other. [@Sun2019b] in another article on the analysis of methylation patterns between CpG sites in a breast cancer study revealed that co-methylation region lengths differed significantly for unmethylated and methylated states. Another finding worthy of note is that by [@sun2022].
They concluded that the co-methylation patterns on chromosome X were different from the other chromosomes suggesting that the degree of co-methylation may differ on different chromosomes. Further more, the aforementioned methods and their respective packages require the user to set many hyperparameters which leads to inconsistent results in a field so important to human life. In 2019, [@Mallik2019a] performed an extensive simulation study on aforementioned methods and their analysis revealed the lack of consistency in the optimal hyperparameter settings. Also, their findings showed that the methods lacked power to detect DMRs. Owing to the recent findings and the shortcomings of the existing methods, we developed the `idDMR`  R package to help the researcher implement the DMR detection method that uses an array-adaptive normalized kernel-weighted statistic developed by [Alhassan et al (D. Alhassan, G.R. Olbricht, A. Adekpedjou, manuscript submitted)]. This method uses a normalized kernel-weighted statistic to borrow co-methylation information from nearby CpG sites to weight site-level statistics. See [Alhassan et al (D. Alhassan, G.R. Olbricht, A. Adekpedjou, manuscript submitted)] for more details. The advantages of our method include: (1) it requires setting only one hyperparameter and hence leads to consistent results, (2)it increases the sensitivity in detecting true DMRs that may exist in less dense regions while still picking up true DMRs in denser regions \par

[idDMR](https://github.com/DanielAlhassan/idDMR) is currently hosted on GitHub. Currently this package contains one DMR detection method in the function `aaDMR()`. This function uses the array-adaptive normalized kernel statistic in (\ref{A}). The primary arguments it accepts include a `CpGsiteAnnotated` object, an adapted object , from the originally created `cpgannotated` object by [@DMRcate2015]; an agglomerate parameter, g, set by the user that is used as the maximum bp distance in collapsing contiguous sites into regions. We describe the key functions in the package and how to use them in the Usage section below.
 
### 
The estimator at the heart of our DMR detection procedure is the array-adaptive normalized kernel-weighted statistic $S(x_{i})$ below due to [Alhassan et al (D. Alhassan, G.R. Olbricht, A. Adekpedjou, manuscript submitted)]:

\begin{align}
S(x_{i}) = Y_{i} + \sum_{j \neq i}^{n}w_{j}(x_{i})Y_{j}
\label{A}
\end{align}
where $Y_{i}$ is the square of `limma`'s t-statistic [@limma] and $w_{j}(x_{i}) = \dfrac{K\left( \dfrac{x_{j}-x_{i}}{h}\right)}{\displaystyle\sum_{j \neq i}^{n}K\left( \dfrac{x_{j}-x_{i}}{h}\right)}$ is a weight measure defined in [Alhassan et al (D. Alhassan, G.R. Olbricht, A. Adekpedjou, manuscript submitted)].

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

# References

