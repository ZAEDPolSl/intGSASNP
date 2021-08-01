# intGSASNP
Integration Gene Set Analysis for SNP.

## Description

intGSASNP is an R package that contains the implementation of six algorithms (ORA, CERNO, MAGENTA, GSEA, iGSEA4GWAS, and GSEA-SNP) selected for studying the importance importance of SNP integration for gene set analysis in GWAS. Features that distinguish intGSASNP include: 
- selection of preferred integration method (minimum p-value method, Fisher's method, Stouffer's method) and permutation method (by Entrez, SNP, or raw SNP),
- possible application of correction for multiple testing,
- option to choose the number of permutations,
- option to apply LD correction.

Additionally, an example of a dataset with sample refSNP ids, entrez ids, and p-values is provided.

## Installation

You can install the package from [GitHub](https://github.com/) with:
``` r
# install.packages("devtools")
devtools::install_github("ZAEDPolSl/intGSASNP")
```

## Manual

Documentation can be accessed using:

``` r
library(intGSASNP)
```
