##==鉴定细胞类型==##
library(SingleR)
BiocManager::install("celldex")
library("celldex")
dir.create("CellType")
refdata <- ImmGenData()
testdata <- GetAssayData(scRNA3, slot="data")
clusters <- scRNA3@meta.data$seurat_clusters
#使用Monaco参考数据库鉴定
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=cellpred$ClusterID, celltype=cellpred$celltype, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_ImmGenData.csv",row.names = F)
scRNA@meta.data$celltype_ImmGen = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_ImmGen'] <- celltype$celltype[i]}
p1 = DimPlot(scRNA, group.by="celltype_ImmGen", repel=T, label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype_ImmGen", repel=T, label=T, label.size=5, reduction='umap')
p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype_ImmGen.png", p1, width=7 ,height=6)
ggsave("CellType/UMAP_celltype_ImmGen.png", p2, width=7 ,height=6)
ggsave("CellType/celltype_ImmGen.png", p3, width=10 ,height=5)
#使用DICE参考数据库鉴定
refdata <- MouseRNAseqData()
# load('~/database/SingleR_ref/ref_DICE_1561s.RData')
# refdata <- ref_DICE
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_MouseRNA.csv",row.names = F)
scRNA@meta.data$celltype_MouseRNA = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_MouseRNA'] <- celltype$celltype[i]}
p4 = DimPlot(scRNA, group.by="celltype_MouseRNA", repel=T, label=T, label.size=5, reduction='tsne')
p5 = DimPlot(scRNA, group.by="celltype_MouseRNA", repel=T, label=T, label.size=5, reduction='umap')
p6 = p3+p4+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype_MouseRNA.png", p4, width=7 ,height=6)
ggsave("CellType/UMAP_celltype_MouseRNA.png", p5, width=7 ,height=6)
ggsave("CellType/celltype_MouseRNA.png", p6, width=10 ,height=5)
#对比两种数据库鉴定的结果
p8 = p1+p4
ggsave("CellType/ImmGen_MouseRNA.png", p8, width=12 ,height=5)

##保存数据
saveRDS(scRNA,'scRNA.rds')

#手动注释
celltype_marker=c(
  "EPCAM",#上皮细胞 epithelial
  "PECAM1",#内皮细胞 endothelial
  "COL3A1",#成纤维细胞 fibroblasts
  "CD163","AIF1",#髓系细胞 myeloid
  "CD79A",#B细胞
  "JCHAIN",#浆细胞 plasma cell
  "CD3D","CD8A","CD4",#T细胞
  "GNLY","NKG7",#NK细胞
  "PTPRC"#免疫细胞
)


celltype_marker=c(
  "EPCAM",#破骨细胞 Osteoclast
  "PECAM1",#内皮细胞 endothelial
  "COL3A1",#成纤维细胞 fibroblasts
  "CD163","AIF1",#髓系细胞 myeloid
  "CD79A",#B细胞
  "JCHAIN",#浆细胞 plasma cell
  "CD3D","CD8A","CD4",#T细胞
  "GNLY","NKG7",#NK细胞
  "PTPRC"#免疫细胞
)

#cellmarker
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
  
  
  "Sost", "Dmp1", "Mepe",#骨细胞 Osteocyte
  "Fgfr3", "Cd44",#干细胞
  "Ctsk", "Mmp9", "Acp5", "Nfatc1",#破骨细胞 Osteoclast
  "Lepr",#Mesenchymal stem cell间充质干细胞
  "Cd105", "Cd166", "Cd44", "Cd73",#Mesenchymal progenitor cell 间充质前体细胞
  "Alpl", "Bglap", "Col1a1", "Ogn", "Runx2", "Sp7"#Osteogenic cell成骨细胞
)

pdf("celltype/marker.featureplot.umap.pdf")
FeaturePlot(scRNA3, features = celltype_marker,reduction="umap")
dev.off()
ggsave("celltype/marker.featureplot.marker.featureplot.umap.png",FeaturePlot(scRNA3, features = celltype_marker,reduction="umap"))
library(readxl)
b=scRNA3@meta.data
marker_celltype=read_xlsx("marker_celltype.xlsx")
pdf("train.marker.dotplot.pdf")
celltype_marker_all=split(marker_celltype$marker,marker_celltype$celltype)
DotPlot(object = scRNA3, features = celltype_marker_all) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
ggsave("train.marker.dotplot.png",DotPlot(object = train, features = unique(top1$gene)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
marker=marker_celltype$marker
DotPlot(object = scRNA3, 
        features=marker, 
        group.by = 'celltype',
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, 
                                                          vjust = 0.5, hjust=0.5))

DotPlot(scRNA3, features = marker_celltype$marker)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
marker_neutrophil="S100a8"
pdf('celltype/neutrophil_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scrna3, 
              features ='S100a8', 
              cols = c("grey", "blue"), 
              reduction = "umap")
dev.off()
ggsave("celltype/neutrophil.umap.png",FeaturePlot(scRNA3, features = marker_neutrophil,reduction="umap"))

marker_myeloid_progenitor="Mpo"
pdf('celltype/myeloid progenitor_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_myeloid_progenitor, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()
ggsave("celltype/myeloid_progenitor.umap.png",FeaturePlot(scRNA3, features = marker_myeloid_progenitor,reduction="umap"))

marker_macrophage="Csf1r"
pdf('macrophage_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_macrophage, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()
ggsave("celltype/macrophage.umap.png",FeaturePlot(scRNA3, features = marker_myeloid_progenitor,reduction="umap"))


marker_dendritic_cell="Siglech"
pdf('dendritic cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_dendritic_cell, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_B_cell="Cd79a"
pdf('B cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_B_cell, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_pro_B_cell="Vpreb1"
pdf('pro_B_cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_pro_B_cell, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_T_cell="Cd3g"
pdf('T_cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_T_cell, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_NK_cell="Klrd1"
pdf('NK_cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_NK_cell, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_NKT_cell=c("Cd3g","Klrd1")
pdf('T_cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_NKT_cell, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_megakaryocyte="Ms4a2"
pdf('megakaryocyte_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_megakaryocyte, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_RED_CELL="Hbb-bt"
pdf('RED_CELL_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA3, 
            features =marker_RED_CELL, 
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_hematopoietic_stem="Cd34"
pdf('hematopoietic_stem_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_hematopoietic_stem, 
            cols = c("grey", "red"), 
            reduction = "tsne")
dev.off()

marker_msc="Col1a1"
pdf('celltype/neutrophil_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scrna3, 
            features ="Col1a1",
            cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

marker_osteocyte=c("Sost", "Dmp1", "Mepe")
pdf('osteocyte_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_osteocyte, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()

marker_stem_cell=c("Fgfr3", "Cd44")
pdf('stem_cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_stem_cell, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()

marker_osteoclast=c("Ctsk", "Mmp9", "Acp5", "Nfatc1")
pdf('osteoclast_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_osteoclast, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()

marker_Mesenchymal_stem_cell="Lepr"
pdf('Mesenchymal_stem_cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_Mesenchymal_stem_cell, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()

marker_Mesenchymal_progenitor_cell=c("Cd105", "Cd166", "Cd44", "Cd73")
pdf('Mesenchymal progenitor cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_Mesenchymal_progenitor_cell, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()

marker_Osteogenic_cell=c("Alpl", "Bglap", "Col1a1", "Ogn", "Runx2", "Sp7")
pdf('Osteogenic cell_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = scRNA, 
            features =marker_Osteogenic_cell, 
            cols = c("grey", "blue"), 
            reduction = "tsne")
dev.off()
VlnPlot(scRNA3,features = celltype_marker,pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = "Ihh",pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = "Siglech",pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = "Hbb-bt",pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = "Col1a1",pt.size = 0,ncol = 1)
VlnPlot(sce_msc,features = "Cd34",pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = "Ms4a2",pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = "Mpo",pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = c("Sost", "Dmp1", "Mepe"),pt.size = 0,ncol = 1)
VlnPlot(scRNA,features = c("Fgfr3", "Cd44"),pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = c("Ctsk", "Mmp9", "Acp5", "Nfatc1"),pt.size = 0,ncol = 1)
VlnPlot(scRNA,features = c("Cd105", "Cd34", "Cd45", "Cd90"),pt.size = 0,ncol = 1)
VlnPlot(scRNA,features = c("Cd105", "Cd166", "Cd44", "Cd73"),pt.size = 0,ncol = 1)
VlnPlot(scRNA3,features = c("Alpl", "Bglap", "Col1a1", "Ogn", "Runx2", "Sp7"),pt.size = 0,ncol = 1)
VlnPlot(scRNA,features = c("", "Bglap", "Col1a1", "Ogn", "Runx2", "Sp7"),pt.size = 0,ncol = 1)
ggsave(filename = "marker.png",device = "png",width = 44,height = 33,units = "cm")

dir.create("enrich")
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

cellpred=read.csv("celltype.csv")
rownames(cellpred)=cellpred$锘緾lusterID
celltype = data.frame(ClusterID=cellpred$锘緾lusterID, celltype=cellpred$celltype, stringsAsFactors = F)
rownames(celltype)=celltype$ClusterID
scRNA3@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA3@meta.data[which(scRNA3@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
library(patchwork)
p1 = DimPlot(scRNA3, group.by="celltype", repel=T, label=T, label.size=3, reduction='tsne')
p2 = DimPlot(scRNA3, group.by="celltype", repel=T, label=T, label.size=3, reduction='umap')
p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype.png", p1, width=7 ,height=6)
ggsave("CellType/UMAP_celltype.png", p2, width=7 ,height=6)
ggsave("CellType/celltype.png", p3, width=10 ,height=5)

sce_msc=subset(scRNA3,idents=c(11,12,13))
sce_macrophage=subset(scRNA3,idents = c(3,4))
saveRDS(sce_macrophage,file = "sce_macrophage.rds")
saveRDS(sce_msc,file = "sce_msc.rds")
saveRDS(scRNA3,file = "scrna3.rds")
