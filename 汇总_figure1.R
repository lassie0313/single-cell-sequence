library(Seurat)
library(tidyverse)
library(patchwork)
dir.create('cluster1')
dir.create('cluster2')
dir.create('cluster3')
set.seed(123) 
dir = c('../DEMO1/filtered_feature_bc_matrix/','../DEMO2/filtered_feature_bc_matrix/','../DEMO3/filtered_feature_bc_matrix/','../DEMO4/filtered_feature_bc_matrix/')
names(dir) = c('S1','S2','S3','S4')
counts <- Read10X(data.dir = dir)
scrna3 = CreateSeuratObject(counts, min.features = 300, min.cells = 3)
dim(scrna3) 
table(scrna3@meta.data$orig.ident)
scrna3 <- NormalizeData(scrna3)
scrna3 <- FindVariableFeatures(scrna3, selection.method = "vst")
scrna3 <- ScaleData(scrna3, features = VariableFeatures(scrna3))
scrna3 <- RunPCA(scrna3, features = VariableFeatures(scrna3))
plot1 <- DimPlot(scrna3, reduction = "pca", group.by="orig.ident")
DimHeatmap(scrna3,
           dims = 1:30,
           cells = 300,
           balanced = TRUE)
plot2 <- ElbowPlot(scrna3, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster1/pca.png", plot = plotc, width = 8, height = 4)
scrna3 <- FindNeighbors(scrna3, dims = pc.num) 
scrna3 <- FindClusters(scrna3, resolution = 0.1)
table(scrna3@meta.data$seurat_clusters)
metadata <- scrna3@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

##非线性降维
#tSNE
scrna3 = RunTSNE(scrna3, dims = pc.num)
embed_tsne <- Embeddings(scrna3, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster1/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scrna3, reduction = "tsne", label=T) 
ggsave("cluster1/tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scrna3, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster1/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster1/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

pc.num=1:30
#UMAP
scrna3 <- RunUMAP(scrna3, dims = pc.num)
embed_umap <- Embeddings(scrna3, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scrna3, reduction = "umap", label=F) 
ggsave("cluster1/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scrna3, reduction = "umap", group.by='orig.ident')
ggsave("cluster1/UMAP.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("cluster1/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)

#合并tSNE与UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)

##==鉴定细胞类型==##
celltype_marker=c(
  "S100a8",#中性粒细胞 neutrophil
  "Mpo",#髓系祖细胞 myeloid progenitor
  "Csf1r",#巨噬细胞 macrophage
  "Siglech",#树突状细胞 dendritic cell
  "Cd79a",#B细胞 B cell
  "Vpreb1",#前B细胞 pro-B cell
  "Cd3g",#T细胞 T cell
  "Klrd1",#NK细胞 natural killer cell
  "Ms4a2",#巨核细胞 megakaryocyte
  "Hbb-bt",#红细胞 erythrocyte
  "Col1a1") #间充质细胞
