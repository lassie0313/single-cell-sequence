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
sce = CreateSeuratObject(counts, min.features = 300, min.cells = 3)
dim(sce) 
table(sce@meta.data$orig.ident)

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.features = 300, min.cells = 3)
}
sce_bone <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], 
                                    scRNAlist[[4]]))
sce_E11_bone <- readRDS("sce_E11_BONE.rds")
sce_E13_cartilage <- readRDS("sce_E13_cartilage.rds")
sce_E13_bone <- readRDS("sce_E13_BONE.rds")
sce_E15_cartilage <- readRDS("sce_E15_cartilage.rds")
sce_E15_bone <- readRDS("sce_E15_BONE.rds")
sce_E18_cartilage <- readRDS("sce_E18_cartilage.rds")
sce_E18_bone <- readRDS("sce_E18_BONE.rds")


sce_bone <- merge(sce_E11_bone, y=c(sce_E13_bone,sce_E15_bone, 
                                    sce_E18_bone), add.cell.ids = c("E11", "E13", "E15", "E18"))

dim(sce_bone)   #查看基因数和细胞总数
table(sce_bone@meta.data$orig.ident)  #查看每个样本的细胞数

scrna3 <- NormalizeData(scrna3)
scrna3 <- FindVariableFeatures(sce, selection.method = "vst")
sce <- ScaleData(sce, features = VariableFeatures(sce))
scrna3 <- RunPCA(scrna3, features = VariableFeatures(scrna3))
plot1 <- DimPlot(sce, reduction = "pca", group.by="orig.ident")
DimHeatmap(sce,
           dims = 1:30,
           cells = 300,
           balanced = TRUE)
plot2 <- ElbowPlot(sce, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster1/pca.png", plot = plotc, width = 8, height = 4)
print(c("请选择哪些pc轴用于后续分析？示例如下：","pc.num=1:15"))
#选取主成分
pc.num=1:23

##细胞聚类
sce_bone <- FindNeighbors(sce_bone, dims = pc.num) 
sce_bone <- FindClusters(sce_bone, resolution = 0.1)
table(sce_bone@meta.data$seurat_clusters)
metadata <- sce_bone@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

##非线性降维
#tSNE
sce_bone = RunTSNE(sce_bone, dims = pc.num)
embed_tsne <- Embeddings(sce_bone, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster1/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(sce_bone, reduction = "tsne", label=T) 
ggsave("cluster1/tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(sce_bone, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster1/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster1/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

pc.num=1:30
#UMAP
sce <- RunUMAP(sce, dims = pc.num)
embed_umap <- Embeddings(sce, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(sce, reduction = "umap", label=F) 
ggsave("cluster1/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(sce, reduction = "umap", group.by='orig.ident')
ggsave("cluster1/UMAP.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("cluster1/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)

#合并tSNE与UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
