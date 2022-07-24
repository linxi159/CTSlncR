# CTSlncR
Inferring cell-type-specific lncRNA regulation in the developing human neocortex with CTSlncR
![](https://github.com/linxi159/CTSlncR/blob/main/figures/Figure_1.tif) 

## Description of each directory and each file
data: The preprocessed data from real scRNA-seq data in GEO.

figures: The plot for CTSlncR.

Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData: Matched lncRNA and mRNA expression data across 276 single cells in the human neocortex, Putative lncRNA-target binding information.

CSlncR_network_bootstrap_GSE71315.RData: 276 cell-specific lncRNA regulatory networks.

CTSlncR.R: Utility functions for exploring cell-type-specific lncRNA regulation.

step1_data_preprocessing_dividing_mRNA_lncRNA.py: Dividing mRNAs and lncRNAs in the scRNA-seq data.

step2_case_study.R: Running scripts for exploring cell-type-specific lncRNA regulation.

## The usage of CTSlncR
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of CTSlncR is implemented in CTSlncR.R. The users can simply run the scripts as follows.
```
* source("step2_case_study.R")
```

## Quick example to use CTSlncR
For identifying cell-type-specific lncRNA regulation, users should prepare lncRNA and mRNA single-cell co-expression data. Paste the datasets and our source file (CTSlncR.R) into a single folder (set the folder as the directory of R environment), users can use the following scripts to identify cell-type-specific lncRNA regulation. For convenience, the datasets prepared for users are from our datasets (Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData).
```
## Load required R packages, please firstly install the following R packages before running scripts
library(pracma)
library(WGCNA)
library(igraph)
library(miRspongeR)
library(biclique)
library(corrplot)
library(dendextend)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(clusterProfiler)
library(msigdbr)
library(vroom)
library(doParallel)

## Load utility functions
setwd("~/CTSlncR")
source("CTSlncR.R")

## Load prepared datasets
ASD <- read.csv("./data/ASD_gene_lncRNAs_mRNAs.csv")
ASD <- as.matrix(ASD)

## Load the preprocessing single-cell sequencing data including filtering genes with expression <= 10 cells, and log e(x+1).
lncRNA_scRNA_raw  <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_lncRNA.csv")
rownames(lncRNA_scRNA_raw ) <- lncRNA_scRNA_raw[,1]
lncRNA_scRNA_raw  <- lncRNA_scRNA_raw[,-1]
lncRNA_scRNA_raw   <- t(lncRNA_scRNA_raw)
lncRNA_scRNA_raw   <- as.data.frame(lncRNA_scRNA_raw)

mRNA_scRNA_raw <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_mRNA.csv")
rownames(mRNA_scRNA_raw) <- mRNA_scRNA_raw[,1]
mRNA_scRNA_raw <- mRNA_scRNA_raw[,-1]
mRNA_scRNA_raw  <- t(mRNA_scRNA_raw)
mRNA_scRNA_raw <- as.data.frame(mRNA_scRNA_raw)

lncRNA_gene_pre <- read.csv("./data/lncRNA_gene-pre.csv")
lncRNA_mRNA_exp <- read.csv("./data/lncRNA_mRNA-exp.csv")
lncRNA_mRNA_exp <- as.matrix(lncRNA_mRNA_exp)

## compute the average expression values of duplicate genes and remove genes with constant expression values in all cells

    # Transformation using log2(x+1)
    lncRNA_scRNA_norm <- log2(lncRNA_scRNA_raw+1)
    mRNA_scRNA_norm <- log2(mRNA_scRNA_raw+1) 

    # Compute the average expression values of duplicate genes
    lncRNA_scRNA_norm_average <- Averg_Duplicate(lncRNA_scRNA_norm)
    mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_scRNA_norm)

    # Remove genes with constant expression values in all cells
    lncRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(lncRNA_scRNA_norm_average)[2]), function(i) sd(lncRNA_scRNA_norm_average[, i])))
    lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_average[, which(lncRNA_scRNA_norm_sd > 0)]
    mRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(mRNA_scRNA_norm_average)[2]), function(i) sd(mRNA_scRNA_norm_average[, i])))
    mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_average[, which(mRNA_scRNA_norm_sd > 0)]

# save    
# save(lncRNA_scRNA_norm, mRNA_scRNA_norm, lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter, lncRNA_gene_pre, lncRNA_mRNA_exp,file = "Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData")
load("Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData")

# # 多核并行计算
# library(parallel)
# #CSlncR_network_bootstrap_null <- CSlncR_net_bootstrap(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, boxsize = 0.1, bootstrap_betw_point = 5, bootstrap_num = 1, p.value.cutoff = 0.01)
# detectCores(logical = F)  # 8
# mc <- getOption("mc.cores", 8)
# system.time({ CSlncR_network_bootstrap_null <- mclapply(list(1), CSlncR_net_bootstrap,lncR=lncRNA_scRNA_norm_filter,mR=mRNA_scRNA_norm_filter, mc.cores = mc);});
# stopCluster(mc);
# 
# #load("CSlncR_network_bootstrap_null_GSE71315.RData")
# CSlncR_network_bootstrap_null = CSlncR_network_bootstrap_null[[1]]

## Discovering cell-specific lncRNA-mRNA regulatory network   
    CSlncR_network_bootstrap_null <- CSlncR_net_bootstrap(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                                          boxsize = 0.1, bootstrap_betw_point = 5, 
                                                          bootstrap_num = 1, p.value.cutoff = 0.01)
```
