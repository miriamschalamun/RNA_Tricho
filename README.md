# RNASeq_analysis
 
## Description
This project contains an R script for analyzing gene expression data using the DESeq2 package for Trichoderma reesei. It was developed as part of Schalamun et al. 2023 The transcription factor STE12 in T. reesei is invovled in carbon and secondary metabolism. The scripts include functions for data preprocessing, normalization, principal component analysis (PCA), heatmaps generation, and gene ontology (GO) enrichment analysis. 
Trichoderma reesei gene annotation is based on PMCID: PMC4771370 and PMC4812632. 
The script was wrritten and executed on Windows 10 and R version 4.2.2. 
 

## Installation

To run the RNASeq_analysis, you will need R installed on your computer. If you do not have R installed, you can download and install it from [CRAN](https://cran.r-project.org/).

Once R is installed, you can run the following commands in your R console to install the required packages:

```R
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "genefilter"))
install.packages(c("readxl", "ggplot2", "dplyr", "ggrepel", "pheatmap", "RColorBrewer", "gplots", "tidyverse", "edgeR", "matrixStats", "xlsx", "dendextend", "topGO", "rrvgo"))
```

## Usage 

1. **Obtain the Script**:
   - Clone this repository to your local machine using the command:
     ```
     git clone https://github.com/yourusername/yourrepository.git
     ```
     Or
   - Download the script file directly from the GitHub repository.

2. **Prepare the Environment**:
   - Launch RStudio or your preferred R environment.
   - Install any required R packages that are not already installed. The script will usually prompt you to install missing packages, or you can find a list of required packages at the beginning of the script.

3. **Set the Working Directory**:
   - Use the `setwd()` function to change your working directory to the location where you have saved the `RNASeq_analysis.Rmd` script. Replace `"/path/to/script"` with the actual file path:
     ```R
     setwd("/path/to/script")
     ```

4. **Run the Script**:
   - Open the `RNASeq_analysis.Rmd` file in RStudio.
   - Run the script interactively by executing code chunks one by one, following the explanations provided within the script. This can be done by clicking the "Run" button within each chunk in RStudio.

5. **Prepare Input Files**:
   - Sample input files specific to *Trichoderma reesei* are provided with this repository. They serve as templates for the format and structure your own data files should have.
   - If you are working with a different organism, you need to create your own input files matching the templates' style and format.
   - Ensure the file names in your script correspond to the names of your actual data files. If you rename your data files, update the file paths and names in the script accordingly.

6. **Review Outputs**:
   - After running the script, check the outputs, which may include processed data files, figures, and logs. Review them to ensure the analysis ran as expected.
  




