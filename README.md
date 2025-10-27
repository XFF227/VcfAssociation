---
title: "VcfAssociation"
output: github_document
---

# VcfAssociation
Title from DESCRIPTION.

## Description
Purpose: Input VCF/phenotype files for association tests between variants/genes, output publishable figures (A3 Q1). Analyzes VCF (merged variant data for individuals), phenotypic CSV/TSV, optional covariates (age/sex; Q2). Improves: Integrates reading/annotation/analysis/visualization, reducing package switching; GWAS module for quick variant ID; model specific variant-phenotype or variant-driver gene associations (Q1). Unique: One-step QC/annotation/single-variant/gene-level tests vs separate vcfR/qqman/VariantAnnotation (Q3). Developed with R 4.3.0 on macOS (sessionInfo()). No Shiny.

## Installation
To install the latest version of the package:
```r
install.packages("devtools")
library("devtools")
devtools::install_github("XFF227/VcfAssociation", build_vignettes = TRUE)
library("VcfAssociation")