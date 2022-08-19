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

This workshop will cover standard bulk RNA-seq analysis using R. For demonstration purposes, we will be using some unpublished prostate cancer data. For time limitations, we will begin will pre-aligned read counts. That is, raw read files (FASTQ files) have already been aligned to the human reference genome to obtain a read counts table. We will load these read counts into R and use DESeq2 to perform differential expression analysis. 

Our experimental data was obtained by performing a drug screen on advanced-stage metastatic prostate cancer biopsies. In brief, tumor biopsies were obtained from five patients and dissociated into cell clumps. Patient tumor cells were seeded into plates and treated with either DMSO or an LSD1 inhibitor (LSD1i or SP2509) for 24 hours. Following this, cells were collected and total RNA was extracted for sequencing on the Illumina HiSeq platform. In addition, clinical attributes from the five patients were collected at the time of biopsy.

You will find read counts in the CSV file ```data/ProstateCancer_DMSO_SP2509_LSD1i_readcounts.csv```. Each column pertains to one sample and each row pertains to a unique transcript. Each patient has two samples: one that is DMSO-treated, and one that is treated with the drug (LDS1i, SP2509).

We also require some sample metadata (also referred to in DESeq2 as ```colData```). This metadata is stored in the CSV file ```data/ProstateCancer_sampleInfo.csv```. This file contains condition information (DMSO or SP2509) and replicate information (replicates 1 to 5 for patients 1 through 5). It also contains some clinical features:

1. mutationCount: total number of mutations detected in the tumor biopsies by whole-genome sequencing
2. diagnosisAge: age of the patient at the time of initial diagnosis (note: not the age at which biopsy was taken)
3. psa: PSA levels at the time of biopsy
4. stage: prostate cancer stage grouping based on the [AJCC (American Joint Committee on Cancer) TNM system](https://cancer.ca/en/cancer-information/cancer-types/prostate/staging).
5. tmb: tumor mutational burden (nonsynonymous mutations only)


## Goals and outline of steps

The goals of this workshop are as follows:

1. Load input data and do some simple data preparation
2. Create DESeq2 data set (dds)
3. Perform differential expression
4. Perform some additional statistical operations
5. Filter the data
6. Make a PCA plot
7. Make a sample clustering plot
8. Make a heatmap to show top differentially expressed genes
9. Make a volcano plot to look at log fold changes and p values
10. Look at some individual gene expression data in relation to clinical attributes

After completing this workshop, you should feel comfortable analyzing bulk RNA-seq data in R!
