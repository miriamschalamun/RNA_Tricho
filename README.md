# RNASeq_analysis
 
## Description
This repository offers an R script for gene expression analysis, tailored for the organism *Trichoderma reesei* and based on the Bioconductor DESeq2 package. For details of the DESeq2 package please refere to the [DESeq2 Vignette]: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html.

The **RNASeq_analysis** script was created to analyze the findings in Schalamun et al. 2023 "The transcription factor STE12 in *T. reesei* is invovled in carbon and secondary metabolism".

The input data is a raw count matrix from *featureCounts* but other raw count matrixes from can be used as well. 
The  scripts include functions for differential gene expression analysis, normalization, principal component analysis (PCA), heatmaps generation, and gene ontology (GO) enrichment analysis. 

*Trichoderma reesei* gene annotation is based on PMCID: PMC4771370 and PMC4812632. 

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

## Example Outputs
I have included some example output files to demonstrate what you can expect from running the RNASeq_analysis script. These examples are based on the *Trichoderma reesei* data and are a good reference for understanding the kind of results the script will produce.

**Count normalization** 
In order to be able to compare counts and visualize them, they need to be normalized first. DESeq2 has its own normalization algorithms, variance stabilizing transformations (VST)  and regularized logarithm (rlog). In this example I use VST normalization. 
Additionally there are other noramlization methods such as counts per million (CPM) and fragments per kilobase (FPKM) that additionally normalizes for fragment lengths. 

 **PCA plot**
A principal component analysis (PCA) plot shows the variation between samples based on their gene expression. Similar samples will cluster together in the plot. The plot has a standard X-Y axis layout, with the axes representing the two principal components that capture the most variation in the data.
   ![PCA Plot Example](output/plots/PCA/example_data_PCA__231115.png)
   Here we see that the main variation of the dataset derives from the different ligth conditions (DD, LL) used (87% variation on x-axis (PC1)). Therefore I split the dataset for the differential gene expression analysis (contrasts) in LL and DD so that the condition does not interfere too much with the effect of mutant vs WT. Here an example PCA for a separated dataset (LL only) and we see that now actually the mutant vs WT are the main variation of the dataset ![PCA Plot Example](output/plots/PCA/example_data_PCA_LL_231115.png)

   The separation of dataset based on the condition is covered in the section "Subset by condition". 

**Heatmaps**
Heatmaps are created using pheatmap and can be used to visualize clustering of samples and genes 
 ![Heatmap Plot Example](output/plots/heatmaps/example_data_heatmap_100_test_231115.png)
Additionally we can focus on certain gene sets as defined in the Anntoation file (T. reesei specific) e.g. only CAZymes or transcription factors. 
![Heatmap Plot Example](output/plots/heatmaps/example_data_heatmap_80__Transcription factors_231115.png)

**Differentially expressed genes** 
DEGs are created using the contrast function which results in two types of file, one in the directory conrasts/all which contains the values for all genes and in the directory contrasts/significant the files are already filtered for p-value and fold change. These values can be changed in the function, normally I used padj < 0.05 and log2 fold change of > |1| (corresonds to a fold change > 2) 

**Annotation**
Gene annotation is performed using the annotation file, this is specific to T. reesei but any other file in the same format can be used. 







  




