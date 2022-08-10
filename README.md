# R data package for IA-SVA 
This R package includes various data examples for [IA-SVA, Lee et. al. (2018)](https://www.nature.com/articles/s41598-018-35365-9). All data simulation code and results used in the manuscript are publicly available at [iasvaExample website](https://dleelab.github.io/iasvaExamples/).


## Installation

To install this package, start R and enter the following commands:

      library(devtools)
      devtools::install_github("dleelab/iasvaExamples")
      

## Load the package
To load this package, enter the following command to the R console:

      library(iasvaExamples)


## Data List


#### Human Islet single cell RNA-seq data from [Lawlor et. al. (2016)](http://genome.cshlp.org/content/early/2017/01/16/gr.212720.116)

Lawlor_Islet_scRNAseq_Read_Counts: Gene-level read counts  
Lawlor_Islet_scRNAseq_Annotations: Annotations describing the islet samples and experimental settings

###### Usage:
      data("Lawlor_Islet_scRNAseq_Read_Counts")
      data("Lawlor_Islet_scRNAseq_Annotations")

#### Human Brain single cell RNA-seq data from [Darmanis et. al. (2015)](http://www.pnas.org/content/112/23/7285.long)

Darmanis_Brain_scRNAseq_Read_Counts: Gene-level read counts  
Darmanis_Brain_scRNAseq_Annotations: Annotations describing the brain samples and experimental settings

###### Usage:
      data("Darmanis_Brain_scRNAseq_Read_Counts")
      data("Darmanis_Brain_scRNAseq_Annotations")

#### Human Islet single cell RNA-seq data from [Xin et. al. (2016)](http://www.cell.com/cell-metabolism/abstract/S1550-4131(16)30434-X)

Xin_Islet_scRNAseq_Read_Counts: Gene-level read counts  
Xin_Islet_scRNAseq_Annotations: Annotations describing the islet samples and experimental settings

###### Usage:
      data("Xin_Islet_scRNAseq_Read_Counts")
      data("Xin_Islet_scRNAseq_Annotations")
