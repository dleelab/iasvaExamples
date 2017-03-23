# data examples for Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) 

## Installation

To install this package, start R and enter the following commands:

      install.packages("devtools")
  
      library(devtools)
  
      install_github("UcarLab/IA-SVA")
      

## RData (.rda) list

### Human Islet single cell RNA-seq data from Lawlor et. al. (2017)

Lawlor_Islet_scRNAseq_Read_Counts: Gene-level read counts  
Lawlor_Islet_scRNAseq_Annotations: Annotations describing the islet samples and experimental settings

#### Usage
data("Lawlor_Islet_scRNAseq_Read_Counts")
data("Lawlor_Islet_scRNAseq_Annotations")
