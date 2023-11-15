# RNASeq_analysis
 
## Description
This project contains R scripts for analyzing gene expression data using the DESeq2 package, among others. It was developed as part of Miriam Schalamun's PhD project on cellular signalling in Trichoderam reesei. The scripts include functions for data preprocessing, normalization, principal component analysis (PCA), heatmaps generation, and gene ontology (GO) enrichment analysis.

## Installation

To run the RNASeq_analysis, you will need R installed on your computer. If you do not have R installed, you can download and install it from [CRAN](https://cran.r-project.org/).

Once R is installed, you can run the following commands in your R console to install the required packages:

```R
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "genefilter"))
install.packages(c("readxl", "ggplot2", "dplyr", "ggrepel", "pheatmap", "RColorBrewer", "gplots", "tidyverse", "edgeR", "matrixStats", "xlsx", "dendextend", "topGO", "rrvgo"))

## Usage

To conduct the gene expression analysis using the provided R scripts, follow the steps below:

1. Clone this repository to your local machine or download the R script files.

2. Open RStudio or another R environment.

3. Set your working directory to the location of the script files. Replace `"/path/to/script"` with the actual path to where you have saved the script:

```R
setwd("/path/to/script")
