## ----load_packages, echo=TRUE--------------------------------------------
rm(list=ls())
library(iasva)
library(iasvaExamples)
library(sva)
library(irlba)
library(Rtsne)
library(pheatmap)
library(corrplot)
library(DescTools) #pcc i.e., Pearson's contingency coefficient
library(RColorBrewer)

color.vec <- brewer.pal(8, "Set1")

#color.pal from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

## ----load_data, echo=TRUE, fig.width=6, fig.height=4---------------------
data("Xin_Islet_scRNAseq_Read_Counts")
data("Xin_Islet_scRNAseq_Annotations")
ls()
counts <- Xin_Islet_scRNAseq_Read_Counts
anns <- Xin_Islet_scRNAseq_Annotations
dim(anns)
dim(counts)

Lib_Size <- colSums(counts)
plot(sort(Lib_Size))
hist(Lib_Size)
summary(Lib_Size)

## ----alpha_cells, echo=TRUE, results='asis'------------------------------
##counts <- counts[, (anns$Cell_Type!="none")] 
##anns <- subset(anns, (Cell_Type!="none"))
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

summary(anns)

Patient_ID <- anns$Donor_ID
Gender <- anns$Gender
Age <- anns$Age
Cell_Type <- anns$Cell_Type
Phenotype <- anns$Condition
Ethnicity <- anns$Ethnicity
Mito_Frac <- anns$Mitochondrial_Fraction

ContCoef(table(Cell_Type, Patient_ID))
ContCoef(table(Cell_Type, Gender))
ContCoef(table(Cell_Type, Age))
ContCoef(table(Cell_Type, Phenotype))
ContCoef(table(Cell_Type, Ethnicity))

## ----geno_lib_size, echo=TRUE, fig.width=6, fig.height=4-----------------
Geo_Lib_Size <- colSums(log(counts+1))
barplot(Geo_Lib_Size, xlab="Cell", las =2)
lcounts <- log(counts + 1)
pca.res = irlba(lcounts - rowMeans(lcounts), 5)$v
cor(Geo_Lib_Size, pca.res[,1])
dim(lcounts)

## ----run_tsne, echo=TRUE, fig.width=7, fig.height=7----------------------
set.seed(34544532)
tsne.res.all <- Rtsne(t(lcounts), dims = 2)

plot(tsne.res.all$Y, main="tSNE", xlab="tSNE Dim1", ylab="tSNE Dim2",pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("topright", levels(Cell_Type), border="white", fill=color.vec, bty="n")

par(mfrow=c(2,2))
plot(tsne.res.all$Y, main="Gender", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Gender], bg=color.vec[Gender], oma=c(4,4,6,12))
legend("topright", levels(Gender), border="white", fill=color.vec, bty="n")
plot(tsne.res.all$Y, main="Patient ID", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=tol21rainbow[Patient_ID], bg=tol21rainbow[Patient_ID], oma=c(4,4,6,12))
legend("topright", levels(Patient_ID), border="white", fill=tol21rainbow, bty="n", cex=0.5)
plot(tsne.res.all$Y, main="Ethnicity", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Ethnicity], bg=color.vec[Ethnicity], oma=c(4,4,6,12))
legend("topright", levels(Ethnicity), border="white", fill=color.vec, bty="n")
plot(tsne.res.all$Y, main="Phenotype", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Phenotype], bg=color.vec[Phenotype], oma=c(4,4,6,12))
legend("topright", levels(Phenotype), border="white", fill=color.vec, bty="n")
par(mfrow=c(1,1))

## ----run_pca, echo=TRUE, fig.width=7, fig.height=6-----------------------
pairs(pca.res[,1:4], main="PCA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], cex=0.8, oma=c(4,4,6,12))
legend("right", levels(Cell_Type), border="white", fill=color.vec, bty="n")

## ----run_sva, echo=TRUE, eval=TRUE, fig.width=7, fig.height=6------------
mod1 <- model.matrix(~Patient_ID+Geo_Lib_Size)
mod0 <- cbind(mod1[,1])
sva.res = svaseq(counts,mod1,mod0, n.sv=4)$sv

pairs(sva.res[,1:4], main="SVA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], cex=0.8, oma=c(4,4,6,12))
legend("right", levels(Cell_Type), border="white", fill=color.vec, bty="n")

## ----run_iasva, echo=TRUE, fig.width=7, fig.height=6---------------------
mod <- model.matrix(~Patient_ID+Geo_Lib_Size)
iasva.res<- iasva(t(counts), mod[,-1],verbose=FALSE, permute=FALSE, num.sv=4)
iasva.sv <- iasva.res$sv

## no color
pairs(iasva.sv[,1:4], pch=21, col="black", bg="black", cex=0.8)

## with color-coding
pairs(iasva.sv[,1:4], main="IA-SVA", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], cex=0.8, oma=c(4,4,6,12))
legend("right", levels(Cell_Type), border="white", fill=color.vec, bty="n")

## ----corr_btw_SVs, echo=TRUE, fig.width=6, fig.height=4------------------
cor(iasva.sv)
corrplot(cor(iasva.sv))

## ----find_markers, echo=TRUE, fig.width=7, fig.height=10-----------------
marker.counts <- find.markers(t(counts), as.matrix(iasva.sv[,c(1,4)]),  rsq.cutoff = 0.3)
nrow(marker.counts)

anno.col <- data.frame(Cell_Type=Cell_Type)
rownames(anno.col) <- colnames(marker.counts)
head(anno.col)
cell.type.col <- color.vec[1:4]
names(cell.type.col) <- c("alpha","beta","delta","PP")
anno.colors <- list(Cell_Type=cell.type.col)

pheatmap(log(marker.counts+1), show_colnames =FALSE, show_rownames = TRUE, clustering_method = "ward.D2",cutree_cols = 5,annotation_col = anno.col, annotation_colors=anno.colors)

## ----run_tsne_post_iasva, echo=TRUE, fig.width=7, fig.height=7-----------
set.seed(75458456)
tsne.res <- Rtsne(unique(t(log(marker.counts+1))), dims = 2)

plot(tsne.res$Y, main="tSNE post IA-SVA", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, col=color.vec[Cell_Type], bg=color.vec[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), border="white", fill=color.vec, bty="n")

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

