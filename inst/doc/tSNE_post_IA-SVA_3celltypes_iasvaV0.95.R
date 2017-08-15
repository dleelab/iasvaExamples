## ----load_packages, echo=TRUE--------------------------------------------
rm(list=ls())
library(irlba)
library(iasva)
library(iasvaExamples)
library(Rtsne)
library(pheatmap)
library(corrplot)
library(DescTools) #pcc i.e., Pearson's contingency coefficient
library(RColorBrewer)

color.vec <- brewer.pal(8, "Set1")

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
counts <- counts[, (anns$Phenotype!="Non-Diabetic")&
                    ((anns$Cell_Type=="GCG")|
                      (anns$Cell_Type=="INS")|
                        (anns$Cell_Type=="KRT19"))]
anns <- subset(anns, (Phenotype!="Non-Diabetic")& 
                      ((Cell_Type=="GCG")|
                        (Cell_Type=="INS")|
                          (Cell_Type=="KRT19")))
dim(counts)
dim(anns)

anns <- droplevels(anns)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

filter = apply(counts, 1, function(x) length(x[x>5])>=3)

counts = counts[filter,]
dim(counts)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

Patient_ID <- anns$Patient_ID
Cell_Type <- anns$Cell_Type
Batch <- anns$Batch

## ----geno_lib_size, echo=TRUE, fig.width=7, fig.height=4-----------------
Geo_Lib_Size <- colSums(log(counts+1))
barplot(Geo_Lib_Size, xlab="Cell", las=2, ylab = "Geometric Library Size")
lcounts <- log(counts + 1)

# PC1 and Geometric library size correlation
pc1 = irlba(lcounts - rowMeans(lcounts), 1)$v[,1] ## partial SVD
cor(Geo_Lib_Size, pc1)

## ----run_tsne, echo=TRUE, fig.width=7, fig.height=7----------------------
set.seed(323542534)
tsne.res <- Rtsne(t(lcounts), dims = 2)

plot(tsne.res$Y, main="tSNE", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), border="white",fill=color.vec, bty="n")

## ----run_iasva, echo=TRUE, fig.width= 7, fig.height=6--------------------
mod <- model.matrix(~Patient_ID+Batch+Geo_Lib_Size)
iasva.res<- iasva(t(counts), mod[,-1],verbose=FALSE, permute=FALSE, num.sv=5)
iasva.sv <- iasva.res$sv

#with color-coding based on true cell-type
pairs(iasva.sv, main="IA-SVA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,14))
legend("right", levels(Cell_Type), border="white", fill=color.vec, bty="n")

cor(Geo_Lib_Size, iasva.sv)

corrplot(cor(iasva.sv))

## ----find_markers, echo=TRUE, fig.width=7, fig.height=10-----------------
marker.counts <- find.markers(t(counts), as.matrix(iasva.sv[,c(1,2)]), rsq.cutoff = 0.4)
nrow(marker.counts)

anno.col <- data.frame(Cell_Type=Cell_Type)
rownames(anno.col) <- colnames(marker.counts)
head(anno.col)

cell.type.col <- color.vec[1:3]
names(cell.type.col) <- c("GCG","INS","KRT19")
anno.colors <- list(Cell_Type=cell.type.col)

pheatmap(log(marker.counts+1), show_colnames =FALSE, fontsize_row = 7, clustering_method = "ward.D2", cutree_cols = 3, annotation_col = anno.col, annotation_colors = anno.colors)

## ----run_tsne_post_iasva, echo=TRUE, fig.width=7, fig.height=7-----------
set.seed(3445462)
tsne.res <- Rtsne(unique(t(log(marker.counts+1))), dims = 2)

plot(tsne.res$Y, main="tSNE post IA-SVA", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), border="white", fill=color.vec, bty="n")

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

