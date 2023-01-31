rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)
scRNA3=readRDS("scrna3.rds")
##第9步：标记基因鉴定和可视化
library(Seurat)
library(ggplot2)
library(dplyr)
markers <- FindAllMarkers(scRNA3,logfc.threshold=0.5,test.use="wilcox",min.pct=0.25,only.pos=TRUE)
head(markers)
##排序，将同一个cluster的marker gene排在一起
markers <- markers %>% group_by(cluster)
write.table(markers,file="cellmarker_pca.xls",sep="\t",row.names=F,col.names=T,quote=F)

## 标记基因可视化
top1 <- markers %>% group_by(cluster) %>% top_n(n = 1, wt=avg_log2FC)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt=avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt=avg_log2FC)
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt=avg_log2FC)
write.table(top10,file="cellmarker_top10.xls",sep="\t",row.names=F,col.names=T,quote=F)
#1, 热图
tmp=scRNA3
scRNA3@assays$RNA@scale.data <- scale(scRNA3@assays$RNA@data, scale = TRUE)
pdf("marker.heatmap.pdf")
DoHeatmap(scRNA3, features = unique(top10$gene),size = 0.5,slot = "scale.data") + NoLegend()

DoHeatmap(scRNA3, features = as.character(unique(top10$gene)), group.by = "seurat_clusters", assay = "RNA", group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B","#FFA500","#24998d", "#a1488e","#3f9337"))+ scale_fill_gradientn(colors = c("navy","white","firebrick3"))  
pdf("marker.heatmap_celltype.pdf",width = 30, height = 10)
DoHeatmap(scRNA3,
          features = as.character(unique(top5$gene)),
          group.by = "celltype",
          assay = 'RNA',
          group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B","#FFA500","#94DA99", "#00BFC4","#3F9337"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
dev.off()
DoHeatmap(scRNA3, features = unique(top10$gene), slot = "data") + NoLegend()
dev.off()
ggsave("marker.heatmap.png",DoHeatmap(scRNA3, features = unique(top2$gene)) + NoLegend())

top1 <- markers %>% group_by(cluster) %>% top_n(n = 1, wt=avg_log2FC)

scRNA3 = RunTSNE(scRNA3, dims = 1:30)
embed_tsne <- Embeddings(scRNA33, 'tsne')


a=c("Itgav", "Cd200", "Cd45", "Thy1", "Cd202b", "Tie2", "Ter119", "6c3", "Cd105")
#2，小提琴图
pdf("marker.vlnplot.pdf")
VlnPlot(scRNA3, features = unique(top1$gene))
dev.off()
ggsave("train.marker.vlnplot.png",VlnPlot(train, features = unique(top1$gene)))

install.packages("remotes")  
remotes::install_github("lyc-1995/MySeuratWrappers")
library(MySeuratWrappers)  
#需要展示的基因  
markers <- c()  
my36colors <-c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87", "#E95C59","#E59CC4","#AB3282", "#23452F","#BD956A","#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD", "#CCE0F5",  "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175")
VlnPlot(scRNA3, features = unique(top1$gene), stacked=T, pt.size=0.9, cols = my36colors, direction = "horizontal", x.lab = "", y.lab = "")+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())



#3，散点图
pdf("train.marker.featureplot.umap.pdf")
FeaturePlot(scRNA3, features = unique(top1$gene),reduction="umap")
dev.off()
ggsave("train.marker.featureplot.umap.png",FeaturePlot(train, features = unique(top1$gene),reduction="umap"))

#4，气泡图
pdf("train.marker.dotplot.pdf")
DotPlot(object = scRNA3, features = unique(top1$gene)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
ggsave("train.marker.dotplot.png",DotPlot(object = train, features = unique(top1$gene)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

##第10步：保存Seurat对象
saveRDS(train, file = "train.rds")
