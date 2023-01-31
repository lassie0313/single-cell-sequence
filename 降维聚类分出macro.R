scRNA3 <- NormalizeData(scRNA3)
scRNA3 <- FindVariableFeatures(scRNA3, selection.method = "vst")
scRNA3 <- ScaleData(scRNA3, features = VariableFeatures(scRNA3))
scRNA3 <- RunPCA(scRNA3, features = VariableFeatures(scRNA3))
plot1 <- DimPlot(scRNA3, reduction = "pca", group.by="orig.ident")
DimHeatmap(scRNA3,
           dims = 1:20,
           cells = 500,
           balanced = TRUE)
plot2 <- ElbowPlot(scRNA3, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster4/pca.png", plot = plotc, width = 8, height = 4)
print(c("请选择哪些pc轴用于后续分析？示例如下：","pc.num=1:15"))
#选取主成分
pc.num=1:30

##细胞聚类
scRNA3 <- FindNeighbors(scRNA3, dims = pc.num) 
scRNA3 <- FindClusters(scRNA3, resolution = 0.02)
table(scRNA3@meta.data$seurat_clusters)
metadata <- scRNA3@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster4/cell_cluster_re0.1.csv',row.names = F)

##非线性降维
#tSNE
scRNA3 = RunTSNE(scRNA3, dims = pc.num)
embed_tsne <- Embeddings(scRNA3, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster4/embed_tsne_re0.1.csv')
#group_by_cluster
plot1 = DimPlot(scRNA3, reduction = "tsne", label=T) 
ggsave("tSNE_re0.15.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNA3, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster4/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster4/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

#UMAP
scRNA3 <- RunUMAP(scRNA3, dims = pc.num)
embed_umap <- Embeddings(scRNA3, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster4/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA3, reduction = "umap", label=T) 
ggsave("UMAP_0.15.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scRNA3, reduction = "umap", group.by='orig.ident')
ggsave("cluster4/UMAP_sample.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("cluster4/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)

#合并tSNE与UMAP
library(patchwork)
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster4/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)

library(Seurat)
library(ggplot2)
library(dplyr)
markers <- FindAllMarkers(scRNA3,logfc.threshold=0.5,test.use="wilcox",min.pct=0.25,only.pos=TRUE)
head(markers)
##排序，将同一个cluster的marker gene排在一起
markers <- markers %>% group_by(cluster)
write.table(markers,file="cluster4/after/cellmarker_pca.xls",sep="\t",row.names=F,col.names=T,quote=F)

## 标记基因可视化
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt=avg_log2FC)
write.table(top10,file="cluster4/after/cellmarker_top10.xls",sep="\t",row.names=F,col.names=T,quote=F)

a=markers[which(markers[,6]==1),]
library(clusterProfiler)
library(org.Mm.eg.db)
# 根据需要更改DEG的值
DEG <- a$gene
#### 第一步，从org.Hs.eg.db提取ENSG的ID 和GI号对应关系
keytypes(org.Mm.eg.db)
degID <- bitr(DEG, fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Mm.eg.db )
head(degID)

enrich <- enrichGO(gene =degID[,2],OrgDb='org.Mm.eg.db',ont="BP",pvalueCutoff=1,qvalueCutoff=1)
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])  ))
enrich_factor <- GeneRatio/BgRatio
out <- data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue,enrich$qvalue,enrich$geneID)
colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue","qvalue","geneID")
# barplot
bar <- barplot(enrich,showCategory=20,title="Biological Pathway",colorBy="p.adjust")
bar
ggsave("enrich/15.png", plot = bar, width = 8, height = 7)

a=c("Ahnak","Mpeg1","Cd68","Csf1r")
VlnPlot(scRNA3,features = a,pt.size = 0,ncol = 1)

marker_macrophage=c("Ahnak","Mpeg1","Cd68","Csf1r")
pdf('macrophage_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_macrophage, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()
ggsave("macrophage.tsne.png",FeaturePlot(scRNA3, features = marker_macrophage,reduction="tsne"))

sce_longbone_macro=subset(scRNA3,idents=c(3))
saveRDS(sce_longbone_macro,file = "sce_longbone_macro.rds")