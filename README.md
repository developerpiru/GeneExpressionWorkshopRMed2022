# GeneExpressionWorkshopRMed2022
Gene expression in R workshop for RMed 2022 conference

## Preamble
This repo contains materials for the Gene Expression Analysis in R workshop for the R/Medicine Conference 2022.

## Setup

Before attending the workshop, please clone this repo and make sure you have a working installation of R with the following packages:
```
# From Cran
install.packages(ggplot2)
install.packages(RColorBrewer)
install.packages(circlize)

# From Bioconductor (requires BiocManager)
BiocManager::install(DESeq2)
BiocManager::install(apeglm)
BiocManager::install(AnnotationDbi)
BiocManager::install(org.Hs.eg.db)

# From GitHub (requires devtools)
devtools::install_github(ComplexHeatmap)
devtools::install_github(EnhancedVolcano)
```

## Introduction

This workshop will cover standard bulk RNA-seq analysis using R. For demonstration purposes, we will be using some unpublished prostate cancer data...

## to do
