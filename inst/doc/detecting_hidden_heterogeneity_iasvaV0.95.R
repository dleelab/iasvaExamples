## ----install_packages, echo=TRUE, eval=FALSE-----------------------------
#  #devtools
#  library(devtools)
#  #iasva
#  devtools::install_github("UcarLab/IA-SVA")
#  #iasvaExamples
#  devtools::install_github("dleelab/iasvaExamples")

## ----load_packages, echo=TRUE, message=FALSE-----------------------------
rm(list=ls())
library(irlba) # partial SVD, the augmented implicitly restarted Lanczos bidiagonalization algorithm
library(iasva)
library(iasvaExamples)
library(sva)
library(Rtsne)
library(pheatmap)
library(corrplot)
library(DescTools) #pcc i.e., Pearson's contingency coefficient
library(RColorBrewer)

color.vec <- brewer.pal(3, "Set1")

## ----load_data, echo=TRUE------------------------------------------------
data("Lawlor_Islet_scRNAseq_Read_Counts")
data("Lawlor_Islet_scRNAseq_Annotations")
ls()
counts <- Lawlor_Islet_scRNAseq_Read_Counts
anns <- Lawlor_Islet_scRNAseq_Annotations
dim(anns)
dim(counts)

summary(anns)
ContCoef(table(anns$Gender, anns$Cell_Type))
ContCoef(table(anns$Phenotype, anns$Cell_Type))
ContCoef(table(anns$Race, anns$Cell_Type))
ContCoef(table(anns$Patient_ID, anns$Cell_Type))
ContCoef(table(anns$Batch, anns$Cell_Type))

## ----alpha_cells, echo=TRUE, results='asis'------------------------------
counts <- counts[, (anns$Phenotype!="Non-Diabetic")&(anns$Cell_Type=="GCG")] 
anns <- subset(anns, (Phenotype!="Non-Diabetic")&(Cell_Type=="GCG"))
dim(counts)
dim(anns)

anns <- droplevels(anns)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

# filter out genes that are sparsely and lowly expressed
filter = apply(counts, 1, function(x) length(x[x>5])>=3)

counts = counts[filter,]
dim(counts)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

## ----geno_lib_size, echo=TRUE, fig.width=7, fig.height=4-----------------
Geo_Lib_Size <- colSums(log(counts+1))
barplot(Geo_Lib_Size, xlab="Cell", ylab="Geometric Lib Size", las=2)
lcounts <- log(counts + 1)

# PC1 and Geometric library size correlation
pc1 = irlba(lcounts - rowMeans(lcounts), 1)$v[,1] ## partial SVD
cor(Geo_Lib_Size, pc1)

## ----run_iasva, echo=TRUE, fig.width= 7, fig.height=6--------------------
set.seed(454353)
Patient_ID <- anns$Patient_ID
mod <- model.matrix(~Patient_ID+Geo_Lib_Size)
iasva.res<- iasva(t(counts), mod[,-1],verbose=FALSE, permute=FALSE, num.sv=5) ##irlba
iasva.sv <- iasva.res$sv

plot(iasva.sv[,1], iasva.sv[,2], xlab="SV1", ylab="SV2")

Cell_Type <- as.factor(iasva.sv[,2] > -0.2) 
levels(Cell_Type)=c("Cell1","Cell2")
table(Cell_Type)

# We identified 6 outlier cells based on SV2 that are marked in red
pairs(iasva.sv, main="IA-SVA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12)) #4,4,6,12
legend("right", levels(Cell_Type), fill=color.vec, bty="n")

plot(iasva.sv[,1:2], main="IA-SVA", pch=21, xlab="SV1", ylab="SV2", col=color.vec[Cell_Type], bg=color.vec[Cell_Type])

cor(Geo_Lib_Size, iasva.sv[,2])

corrplot(cor(iasva.sv))

## ----find_markers, echo=TRUE, fig.width=7, fig.height=14-----------------
marker.counts <- find.markers(t(counts), as.matrix(iasva.sv[,2]))
nrow(marker.counts)
rownames(marker.counts)

anno.col <- data.frame(Cell_Type=Cell_Type)
rownames(anno.col) <- colnames(marker.counts)
head(anno.col)

pheatmap(log(marker.counts+1), show_colnames =FALSE, clustering_method = "ward.D2",cutree_cols = 2,annotation_col = anno.col)

## ----run_tsne, echo=TRUE, fig.width=7, fig.height=7----------------------
set.seed(323542534)
tsne.res <- Rtsne(t(lcounts), dims = 2)

plot(tsne.res$Y, main="tSNE", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), fill=color.vec, bty="n")

## ----run_tsne_post_iasva, echo=TRUE, fig.width=7, fig.height=7-----------
set.seed(345233)
tsne.res <- Rtsne(unique(t(log(marker.counts+1))), dims = 2)

plot(tsne.res$Y, main="tSNE post IA-SVA", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), fill=color.vec, bty="n")

## ----run_pca, echo=TRUE, fig.width=7, fig.height=6-----------------------
pca.res = irlba(lcounts - rowMeans(lcounts), 5)$v ## partial SVD

pairs(pca.res, main="PCA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12)) #4,4,6,12
legend("right", levels(Cell_Type), fill=color.vec, bty="n")

plot(pca.res[,2:3], main="PCA", xlab="PC2", ylab="PC3", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), fill=color.vec, bty="n")

## ----run_sva, echo=TRUE, fig.width=7, fig.height=6-----------------------
mod1 <- model.matrix(~Patient_ID+Geo_Lib_Size)
mod0 <- cbind(mod1[,1])

sva.res = svaseq(counts,mod1,mod0, n.sv=5)$sv

pairs(sva.res, main="SVA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12)) #4,4,6,12
legend("right", levels(Cell_Type), fill=color.vec, bty="n")

plot(sva.res[,1:2], main="SVA", xlab="SV1", ylab="SV2", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type])

## ----SV2_geometric_lib_size, echo=TRUE-----------------------------------
cor(Geo_Lib_Size, iasva.sv[,2])

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

