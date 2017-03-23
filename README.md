# Data examples for Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) 





## Installation

To install this package, start R and enter the following commands:

      if(!require(devtools)){
            install.packages("devtools")
            library(devtools)
      }
      install_github("dleelab/iasvaExamples")
      




## Data List


#### Human Islet single cell RNA-seq data from [Lawlor et. al. (2016)](http://genome.cshlp.org/content/early/2017/01/16/gr.212720.116)

Lawlor_Islet_scRNAseq_Read_Counts: Gene-level read counts  
Lawlor_Islet_scRNAseq_Annotations: Annotations describing the islet samples and experimental settings

###### Usage:
      data("Lawlor_Islet_scRNAseq_Read_Counts")
      data("Lawlor_Islet_scRNAseq_Annotations")
