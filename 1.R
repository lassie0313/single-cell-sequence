rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Seurat)
#设置参数：字符型不默认转变为因子型
options(stringAsFactors=F)
#读入cellranger的定量结果
data=Read10X("filtered_feature_bc_matrix/")
sce_S1 <- CreateSeuratObject(data,"S1")
sce_S1

sce_S1@meta.data$group <- "S1"
head(colnames(sce_S1@assays$RNA)) #也就是视频里说的sce@data

sce_S1 <- NormalizeData(sce_S1)
sce_S1 <- ScaleData(sce_S1, display.progress = F)

sce_S1<- FindVariableFeatures(sce_S1, do.plot = F)

a=VariableFeatures(sce_S1)
g.1 <- head(VariableFeatures(sce_S1),1000)
g.2 <- head(VariableFeatures(sce_S2),1000)
g.3 <- head(VariableFeatures(sce_S3),1000)
g.4 <- head(VariableFeatures(sce_S4),1000)

gene.use <- unique(c(g.1 ,g.2 ,g.3 ,g.4))
a=head(rownames(sce_E11@assays$RNA@scale.data))
gene.use <- intersect(gene.use, rownames(sce_S1@assays$RNA@scale.data))
gene.use <- intersect(gene.use, rownames(sce_S2@assays$RNA@scale.data))
gene.use <- intersect(gene.use, rownames(sce_S3@assays$RNA@scale.data))                     
gene.use <- intersect(gene.use, rownames(sce_S4@assays$RNA@scale.data))                      
head(gene.use)
length(gene.use)

colnames(sce_S1@assays$RNA@scale.data) <-  paste0("S1.",colnames(sce_S1@assays$RNA))

start_time <- Sys.time()
sce.comb <- RunMultiCCA(list(sce_S1,sce_S2,sce_S3,sce_S4),
                   add.cell.ids=c("S1.","S2.","S3.","S4."),
                   gene.use=gene.use,
                   num.cc = 30)
a <- RunCCA(object1 = sce_E11,object2 = sce_E13)
