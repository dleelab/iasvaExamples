rm(list=ls())
library(SummarizedExperiment)

pwd <- "/Users/leed1/Desktop/My_Project/GitHub_Projects/iasvaExamples/"
pwd.data <- "/Users/leed1/Desktop/My_Project/Data/Darmanis_Brain_recount2_SRP057196/"
pwd.out <- paste0(pwd,"data/")

gene.table <- read.table("/Users/leed1/Desktop/My_Project/Data/gene_table_ensg/HG19_v74_FeatureData.csv", header=TRUE, sep=",")
dim(gene.table)
table(gene.table$Gene.Biotype)
gene.table <- subset(gene.table, Gene.Biotype=="protein_coding"|Gene.Biotype=="lincRNA")
dim(gene.table)

load(paste0(pwd.data,"rse_gene.Rdata"))
str(rse_gene)
class(rse_gene)
colData(rse_gene)

## make counts
counts <- assay(rse_gene, "counts")
rod<- rowData(rse_gene)
rownames(counts) <- names(rod$symbol)
counts <- counts[rownames(counts)%in%gene.table$X, ] ## use only protein coding genes and lincRNAs.
dim(counts) #nrow: 25537, ncol:461
dropout.rate <- sum(counts==0)/(dim(counts)[1]*dim(counts)[2])*100
dropout.rate
filter = apply(counts, 1, function(x) length(x[x>5])>=3)
counts = counts[filter,]
dropout.rate <- sum(counts==0)/(dim(counts)[1]*dim(counts)[2])*100
dropout.rate
dim(counts) # nrow: 23326, ncol:461

## make meta data
meta <- read.table(paste0(pwd.data, "SRP057196.tsv"), header=TRUE, sep="\t")
dim(meta)   # nrow: 466, ncol:21
meta <- meta[meta$run%in%as.factor(colnames(counts)),]
dim(counts) # nrow: 25537, ncol:461
dim(meta)   # nrow: 461, ncol:21
# extract meta data
meta.list <- lapply(as.character(meta$characteristics), function(x) strsplit(x, ","))
meta.df <- data.frame(matrix(unlist(meta.list), nrow=nrow(meta), byrow=T))
colnames(meta.df) <- c("tissue","cell_type","age","c1_chip_id", "sample_name")
meta.df$run <- meta$run
meta.df$tissue <- as.factor(gsub(".*: ", "", meta.df$tissue))
meta.df$cell_type <- as.factor(gsub(".*: ", "", meta.df$cell_type))
meta.df$age <- as.factor(gsub(".*: ", "", meta.df$age))
meta.df$c1_chip_id <- as.factor(gsub(".*: ", "", meta.df$c1_chip_id))
meta.df$sample_name <- as.factor(gsub(".*: |)", "", meta.df$sample_name))
summary(meta.df)
meta.df <- meta.df[,c(6,1,2,3,4,5)]

## prep data
Darmanis_Brain_scRNAseq_Read_Counts  <- counts
Darmanis_Brain_scRNAseq_Annotations <- meta.df
save(Darmanis_Brain_scRNAseq_Read_Counts, file=paste0(pwd.out,"Darmanis_Brain_scRNAseq_Read_Counts.rda"))
save(Darmanis_Brain_scRNAseq_Annotations, file=paste0(pwd.out,"Darmanis_Brain_scRNAseq_Annotations.rda"))

