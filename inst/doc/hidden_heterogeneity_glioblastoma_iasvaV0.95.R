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
data("Patel_Glioblastoma_scRNAseq_Read_Counts")
data("Patel_Glioblastoma_scRNAseq_Annotations")
ls()
counts <- Patel_Glioblastoma_scRNAseq_Read_Counts
anns <- Patel_Glioblastoma_scRNAseq_Annotations
dim(anns)
dim(counts)

summary(anns)
table(anns$patient_id, anns$subtype)
ContCoef(table(anns$patient_id, anns$subtype))

## ----MGH30_cells, echo=TRUE, results='asis'------------------------------
counts <- counts[, (anns$subtype!="None")&(anns$patient_id=="MGH30")] 
anns <- subset(anns, (subtype!="None")&(patient_id=="MGH30"))
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


Subtype <- anns$subtype
Patient_ID <- anns$patient_id

## ----geno_lib_size, echo=TRUE, fig.width=7, fig.height=4-----------------
Geo_Lib_Size <- colSums(log(counts+1))
barplot(Geo_Lib_Size, xlab="Cell", ylab="Geometric Lib Size", las=2)
lcounts <- log(counts + 1)
pc1 = irlba(lcounts - rowMeans(lcounts), 1)$v[,1] ## partial SVD
cor(Geo_Lib_Size, pc1)

## ----run_iasva, echo=TRUE, fig.width= 7, fig.height=6--------------------
set.seed(345)
mod <- model.matrix(~Geo_Lib_Size)
iasva.res<- iasva(t(counts), mod[,-1],verbose=FALSE, permute=FALSE, num.sv=5) ## irlba
iasva.sv <- iasva.res$sv

Cell_Cycle <- as.factor(iasva.sv[,2] > -0.1) 
levels(Cell_Cycle)=c("Cycle1","Cycle2")
table(Cell_Cycle)

pairs(iasva.sv[,1:5], main="IA-SVA", pch=21, col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle], oma=c(4,4,6,14))
legend("right", levels(Cell_Cycle), fill=color.vec, bty="n")

plot(iasva.sv[,1:2], main="IA-SVA", pch=21, xlab="SV1", ylab="SV2", col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle])
#legend("bottomright", levels(Cell_Cycle), fill=color.vec, bty="n")

cor(Geo_Lib_Size, iasva.sv[,2])

corrplot(cor(iasva.sv))

## ----find_markers, echo=TRUE, fig.width=7, fig.height=12-----------------
marker.counts <- find.markers(t(counts), as.matrix(iasva.sv[,2]), rsq.cutoff = 0.4)
nrow(marker.counts) #87 58
rownames(marker.counts)
anno.col <- data.frame(Subtype=Subtype, Cell_Cycle=Cell_Cycle, Lib_Size=colSums(counts))
rownames(anno.col) <- colnames(marker.counts)
head(anno.col)

pheatmap(log(marker.counts+1), show_colnames =FALSE, clustering_method = "ward.D2",cutree_cols = 2,annotation_col = anno.col)

## ----run_tsne, echo=TRUE, fig.width=7, fig.height=7----------------------
set.seed(323542534)
tsne.res <- Rtsne(t(lcounts), dims = 2, perplexity = 15)

plot(tsne.res$Y, main="tSNE", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Cycle), fill=color.vec, bty="n")

## ----run_tsne_post_iasva, echo=TRUE, fig.width=7, fig.height=7-----------
set.seed(34523)
tsne.res <- Rtsne(unique(t(log(marker.counts+1))), dims = 2, perplexity = 15)

plot(tsne.res$Y, main="tSNE post IA-SVA", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Cycle), fill=color.vec, bty="n")

## ----pca_plot, echo=TRUE, fig.width=7, fig.height=6----------------------
pca.res = irlba(lcounts - rowMeans(lcounts), 5)$v ## partial SVD

pairs(pca.res[,1:5], main="PCA", pch=21, col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle],
      oma=c(4,4,6,14))
legend("right", levels(Cell_Cycle), fill=color.vec, bty="n")

plot(pca.res[,1:2], main="PCA", pch=21, xlab="PC1", ylab="PC2", col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle])

## ----run_sva, echo=TRUE, fig.width=7, fig.height=6-----------------------
mod1 <- model.matrix(~Geo_Lib_Size)
mod0 <- cbind(mod1[,1])

sva.res = svaseq(counts,mod1,mod0, n.sv=5)$sv

pairs(sva.res[,1:5], main="SVA", pch=21, col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle], oma=c(4,4,6,12)) #4,4,6,12
legend("right", levels(Cell_Cycle), fill=color.vec, bty="n")

plot(sva.res[,1:2], main="SVA", xlab="SV1", ylab="SV2", pch=21, col=color.vec[Cell_Cycle], bg=color.vec[Cell_Cycle])

## ----SV2_geometric_lib_size, echo=TRUE-----------------------------------
cor(Geo_Lib_Size, iasva.sv[,2])

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

