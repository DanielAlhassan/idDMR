---
title: 'idDMR: Identify Differentially Methylated Regions for Microarray Data '
tags:
  - R
  - bioinformatics
  - differentially methylated region
  - kernel smoothing
  - Illumina
authors:
  - name: Daniel Alhassan
    corresponding: true 
    orcid: 0000-0003-0944-5818
    equal-contrib: true
    affiliation: 1 
  - name: Gayla Olbricht
    orcid: 0000-0002-1213-2241
    equal-contrib: true 
    affiliation: 1
  - name: Akim Adekpedjou
    orcid: 0000-0001-9584-4297
    equal-contrib: true 
    affiliation: 1
  - name: Ebenezer Agbozo
    orcid: 0000-0002-2413-3815
    equal-contrib: true 
    affiliation: 2
affiliations:
 - name: Department of Mathematics and Statistics, Missouri University of Science and Technology, USA
   index: 1
 - name: Department of Big Data Analytics and Methods of VideoAnalytics, Institute of               Radioelectronics and Information Technology, Ural Federal University
   index: 2
date: 06 June 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
`idDMR` is an R package currently hosted in GitHub repository for identifying or detecting differentially methylated region in the Human genome for data  [idDMR](https://github.com/DanielAlhassan/idDMR). 


# Statement of need
 Identifying differentially methylated regions is a common but challenging task the researcher faces owing to the complex nature and degree of co-methylation between neighboring CpG sites on a chromosome. Several efforts have been made, methods have been proposed and software packages written for these methods to make the users task easier. As our knowledge on the degree of co-methylation increases, better methods that utilize the new findings are proposed and hence software need to be created to make them available to the user to implement. 

Some existing R packages in this area of research are the  `bumphunter` pacakge [@bumphunter;@minfi], the `ChAMP` package with implements the `Probe Lasso` method to DMR identification
 
 As with the above R/Biconductor packages,  `idDMR` package was designed to allow the user to implement the methods discussed in [Alhassan et 2021]. At the heart of the detection procedure is the use of the normalized kernel-weighted statistic to incorporate co-methylation information at the CpG site-level. Our modeling framework borrows it's foundation work from DMRcate [@DMRcate2015]. 
 

# Mathematics
The estimator at the heart of our DMR detection procedure is the array-adaptive normalized kernel-weighted statistic $S(x_{i})$ below:

$$\begin{align}
S(x_{i})   & = Y_{i} + \sum_{\substack{j = 1\\j \neq i}}^{n}w_{j}(x_{i})Y_{j} 
\end{align}$$

where $Y_{i}$ is the square of `limma`'s t-statistic [@limma] and $w_{j} = \dfrac{K\left( \dfrac{x_{j}-x_{i}}{h}\right)}{\displaystyle\sum_{\substack{j =1\\j \neq i}}^{n}K\left( \dfrac{x_{j}-x_{i}}{h}\right)}$. 
See [Alhassan] for a more in dept discussion of the methods.




# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for .

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"



# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References