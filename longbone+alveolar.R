library(Seurat)
sce_longbone_macro=readRDS("sce_longbone_macro.rds")
scRNA_alver_macro=readRDS("../sce_macrophage.rds")
sce_longbone_macro[["percent.mt"]] <- PercentageFeatureSet(sce_longbone_macro, pattern = "^Mt")

#计算ercc
sce_longbone_macro[["percent.ercc"]] <- PercentageFeatureSet(sce_longbone_macro, pattern = "^Ercc")
sce_longbone_macro[["percent.rp"]] <- PercentageFeatureSet(sce_longbone_macro, pattern = "^Rp")
a=sce_longbone_macro@meta.data
all.genes <- rownames(sce_longbone_macro)
sce_longbone_macro <- ScaleData(sce_longbone_macro, features = all.genes,
                 vars.to.regress = c("percent.rp","percent.mt"))

library(ggplot2)
minGene=500
pctMT=0.5

##数据质控
sce_longbone_macro <- subset(sce_longbone_macro, subset = nFeature_RNA > minGene)
col.num <- length(levels(as.factor(sce_longbone_macro@meta.data$orig.ident)))
violin <-VlnPlot(sce_longbone_macro, group.by = "orig.ident",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc","percent.rp"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 



scRNA_macro_two <- merge(scRNA_alver_macro, y=sce_longbone_macro,add.cell.ids = c("alveolar", "longbone"))


cellpred=read.csv("position.csv")
position = data.frame(orig.ident=cellpred$锘縪rig.ident, position=cellpred$position, stringsAsFactors = F)
scRNA_macro_two@meta.data$position="NA"

for(i in 1:nrow(position)){
  scRNA_macro_two@meta.data[which(scRNA_macro_two@meta.data$orig.ident == position$orig.ident[i]),'position'] <- position$position[i]}
scRNA_macro_two@meta.data$orig.ident=scRNA_macro_two@meta.data$position
a=scRNA_macro_two@meta.data


scRNAlist <- SplitObject(scRNA_macro_two,split.by = "orig.ident")
for(i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]],verbose = FALSE)
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = "vst",nfeatures = 600,
                                         verbose = FALSE)
}
reference.list <- scRNAlist[c("alveolar","longbone")]
scRNA.anchors <- FindIntegrationAnchors(object.list = reference.list,dims = 1:30)
scrna.integrated <- IntegrateData(anchorset = scRNA.anchors,dims = 1:30)
#切换到整合后的assay
DefaultAssay(scrna.integrated) <- "integrated"
scRNA_macro_two=scrna.integrated

#运行标准化可视化与聚类
library(Seurat)
scRNA_macro_two <- NormalizeData(scRNA_macro_two)
scRNA_macro_two <- FindVariableFeatures(scRNA_macro_two, selection.method = "vst")
scRNA_macro_two <- ScaleData(scRNA_macro_two, features = VariableFeatures(scRNA_macro_two))
scRNA_macro_two <- RunPCA(scRNA_macro_two, features = VariableFeatures(scRNA_macro_two))

scRNA_macro_two = RunTSNE(scRNA_macro_two, dims = 1:30)
scRNA_macro_two <- RunUMAP(scRNA_macro_two, reduction = "pca", dims = 1:30)
scRNA_macro_two <- FindNeighbors(scRNA_macro_two, reduction = "pca", dims = 1:30)
scRNA_macro_two <- FindClusters(scRNA_macro_two, resolution = 0.05)

library(ggplot2)
#可视化
p1 <- DimPlot(scRNA_macro_two, reduction = "umap", group.by = "position")
p2 <- DimPlot(scRNA_macro_two, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
p3=DimPlot(scRNA_macro_two, reduction = "umap", split.by = "position")
ggsave("two/umap_split.png", plot = p3, width = 8, height = 4)
ggsave("two/umap_umap.png", plot = p2, width = 8, height = 4)

p1 <- DimPlot(scRNA_macro_two, reduction = "tsne", group.by = "position")
p2 <- DimPlot(scRNA_macro_two, reduction = "tsne", label = TRUE, repel = TRUE)
p1 + p2
p3=DimPlot(scRNA_macro_two, reduction = "tsne", split.by = "position")
ggsave("two/tsne_split.png", plot = p3, width = 8, height = 4)
ggsave("two/tsne_umap.png", plot = p2, width = 8, height = 4)

scRNA_macro_two <- FindNeighbors(scRNA_macro_two, dims = pc.num) 
scRNA_macro_two <- FindClusters(scRNA_macro_two, resolution = 0.001)
table(scRNA_macro_two@meta.data$seurat_clusters)
metadata <- scRNA_macro_two@meta.data
#鉴定不同刺激处理下的差异表达基因
DefaultAssay(scRNA_macro_two) <- "RNA"
#需要先安装包
BiocManager::install('multtest')
install.packages('metap')

nk.markers <- FindConservedMarkers(scRNA_macro_two, grouping.var = "position", verbose = FALSE)
head(nk.markers)

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

scRNA_macro_two <- RenameIdents(scRNA_macro_two, `0` = "alveolar", `1` = "alveolar", `2` = "alveolar", `3` = "alveolar", `4` = "longbone", `5` = "longbone", `6` = "longbone")
table(Idents(scRNA_macro_two))

Idents(scRNA_macro_two) <- "position"
avg <- as.data.frame(log1p(AverageExpression(scRNA_macro_two, verbose = FALSE)$RNA))
avg$gene <- rownames(avg)

genes.to.label = c("Osm", "Grn",  "Tgfb1", "Vegfa", "Adgre5", "Spn", "Cd83", "Mif", "Sema4d", "Cklf", "Copa", "Nrg1")																						

a=c("Osm", "Tgfb1", "Vegfa", "Cd83", "Nrg1")
p1 <- ggplot(avg, aes(alveolar, longbone)) + geom_point() + ggtitle("macrophage")
p1 <- LabelPoints(plot = p1, points = a, repel = TRUE)
ggsave("two/marker.png", plot = p1, width = 8, height = 4)
#使用FindMarkers函数寻找差异表达基因 
library(patchwork)
DEG <- FindMarkers(scRNA_macro_two, ident.1 = "alveolar", ident.2 = "longbone", verbose = FALSE)
DEG$gene=rownames(DEG)
head(b.interferon.response, n = 15)

genes.to.label = c("Sema4d", "Cklf", "Copa", "Nrg1")

p3=FeaturePlot(scRNA_macro_two, features ="Itgal", split.by = "position", max.cutoff = 3, 
            cols = c("grey", "red"))

plots <- VlnPlot(scRNA_macro_two, features = "Itgal", split.by = "position", group.by = "position", 
                 pt.size = 0, combine = FALSE)
plots=wrap_plots(plots = plots, ncol = 1)
ggsave("two/1.png", plot = p1ots, width = 8, height = 4)
saveRDS(scRNA_macro_two,file = "scRNA_macro_two.rds")

#final
a=c("Osm", "Tgfb1", "Vegfa", "Cd83", "Nrg1")
p3=FeaturePlot(scRNA_macro_two, features =a, split.by = "position", max.cutoff = 3, 
               cols = c("grey", "red"))
plots <- VlnPlot(scRNA_macro_two, features = a, split.by = "position", group.by = "position", 
                 pt.size = 0, combine = FALSE)
plots=wrap_plots(plots = plots, ncol = 1)


saveRDS(scRNA_macro_two,file="scRNA_macro_two.rds")
