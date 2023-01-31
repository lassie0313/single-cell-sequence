rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
sce_msc <- readRDS("sce_msc.rds")
library(ggplot2)
sce_msc <- NormalizeData(sce_msc)
sce_msc <- FindVariableFeatures(sce_msc, selection.method = "vst")
sce_msc <- ScaleData(sce_msc, features = VariableFeatures(sce_msc))
sce_msc <- RunPCA(sce_msc, features = VariableFeatures(sce_msc))
plot1 <- DimPlot(sce_msc, reduction = "pca", group.by="orig.ident")
DimHeatmap(sce_msc,
           dims = 1:10,
           cells = 500,
           balanced = TRUE)
plot2 <- ElbowPlot(sce_msc, ndims=10, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster4/pca.png", plot = plotc, width = 8, height = 4)
print(c("请选择哪些pc轴用于后续分析？示例如下：","pc.num=1:15"))
#选取主成分
pc.num=1:30

##细胞聚类
sce_msc <- FindNeighbors(sce_msc, dims = pc.num) 
sce_msc <- FindClusters(sce_msc, resolution = 0.1)
table(sce_msc@meta.data$seurat_clusters)
metadata <- sce_msc@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster4/cell_cluster_re0.1.csv',row.names = F)

##非线性降维
#tSNE
sce_msc = RunTSNE(sce_msc, dims = pc.num)
embed_tsne <- Embeddings(sce_msc, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster4/embed_tsne_re0.1.csv')
#group_by_cluster
plot1 = DimPlot(sce_msc, reduction = "tsne", label=T) 
ggsave("cluster4/tSNE_re0.1.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(sce_msc, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster4/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster4/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

#UMAP
sce_msc <- RunUMAP(sce_msc, dims = pc.num)
embed_umap <- Embeddings(sce_msc, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster4/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(sce_msc, reduction = "umap", label=T) 
ggsave("cluster4/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(sce_msc, reduction = "umap", group.by='orig.ident')
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
markers <- FindAllMarkers(sce_msc,logfc.threshold=0.5,test.use="wilcox",min.pct=0.25,only.pos=TRUE)
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
bar2 <- barplot(enrich,showCategory=20,title="Biological Pathway",colorBy="p.adjust")
bar2
#1, 热图

sce_msc@assays$RNA@scale.data <- scale(sce_msc@assays$RNA@data, scale = TRUE)
pdf("cluster4/after/marker.heatmap.pdf")
DoHeatmap(sce_msc, features = unique(top2$gene),size = 0.5,slot = "scale.data") + NoLegend()

top1 <- markers %>% group_by(cluster) %>% top_n(n = 1, wt=avg_log2FC)

a=c("Bglap", "Lepr", "Cdh5", "Plp1")
#2，小提琴图
pdf("cluster4/after/marker.vlnplot.pdf")
VlnPlot(sce_msc, features = a,ncol = 4)
dev.off()
ggsave("cluster4/after/marker.vlnplot.png",VlnPlot(sce_msc, features = unique(top1$gene)))

#3，散点图
pdf("cluster4/after/marker.featureplot.umap.pdf")
FeaturePlot(sce_msc, features = a,reduction="umap")
dev.off()
ggsave("msc/marker.featureplot.umap.png",FeaturePlot(sce_msc, features = a,reduction="umap"))

#4，气泡图
pdf("cluster4/after/marker.dotplot.pdf")
DotPlot(object = sce_msc, features = a) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
ggsave("msc/marker.dotplot.png",DotPlot(object = sce_msc, features = a) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

##第10步：保存Seurat对象
saveRDS(sce_msc, file = "sce_msc.rds")

a=markers[which(markers[,6]==3),]
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
ggsave("cluster4/after/8.png", plot = bar, width = 8, height = 7)

VlnPlot(sce_msc,features = "Cd206",pt.size = 0,ncol = 1)
VlnPlot(sce_msc,features = "Cd86",pt.size = 0,ncol = 1)
VlnPlot(sce_msc,features = "Cd80",pt.size = 0,ncol = 1)
VlnPlot(sce_msc,features = "Plp1",pt.size = 0,ncol = 1)
a=c("Lepr","Bglap","Cdh5","Plp1")
pdf('B cell_FeaturePlot.pdf', width=10, height=15)
p1=FeaturePlot(object = sce_msc, 
               features ="Plp1", 
               cols = c("grey", "blue"), 
               reduction = "umap")
ggsave("msc/Plp1.png",plot = p1, width=8,height = 7)

cellpred=read_xlsx("msccelltype.xlsx")
rownames(cellpred)=cellpred$ClusterID
celltype = data.frame(ClusterID=cellpred$ClusterID, celltype=cellpred$celltype, stringsAsFactors = F)
rownames(celltype)=celltype$ClusterID
sce_msc@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce_msc@meta.data[which(sce_msc@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
a=sce_msc@meta.data
library(patchwork)
p1 = DimPlot(sce_msc, group.by="celltype", repel=T, label=T, label.size=3, reduction='tsne')
p2 = DimPlot(sce_msc, group.by="celltype", repel=T, label=T, label.size=3, reduction='umap')
p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("msc/tSNE_celltype2.png", p1, width=7 ,height=6)
ggsave("msc/UMAP_celltype2.png", p2, width=7 ,height=6)
ggsave("msc/celltype2.png", p3, width=10 ,height=5)

DotPlot(object = sce_msc, features = a, group.by = 'celltype' ,
        + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
DotPlot(object = sce_msc, 
        features=a, 
        group.by = 'celltype',
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, 
                                                          vjust = 0.5, hjust=0.5))
