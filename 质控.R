##==数据质控==#
scRNA3 <- scRNA3  #以后的分析使用整合的数据进行
##meta.data添加信息
proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA3)))
rownames(proj_name) <- row.names(scRNA3@meta.data)
scRNA3 <- AddMetaData(scRNA3, proj_name)

##切换数据集
DefaultAssay(scRNA3) <- "RNA"
dir.create('QC')
##计算线粒体和红细胞基因比例
scRNA3[["percent.mt"]] <- PercentageFeatureSet(scRNA3, pattern = "^Mt")

#计算ercc
scRNA3[["percent.ercc"]] <- PercentageFeatureSet(scRNA3, pattern = "^Ercc")
library(ggplot2)
violin=VlnPlot(object = scRNA3, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.ercc","percent.mt"), 
        group.by = 'orig.ident',
        ncol = 4)
ggsave("QC/vlnplot_group_before_qc.pdf", plot = violin, width = 12, height = 6)
ggsave("QC/vlnplot_group_before_qc.png", plot = violin, width = 12, height = 6) 
#head(scRNA3@meta.data)
col.num <- length(levels(as.factor(scRNA3@meta.data$orig.ident)))

##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
violin <-VlnPlot(scRNA3, group.by = "proj_name",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA3, feature1 = "nCount_RNA", feature2 = "percent.ercc")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

##设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=500
pctMT=0.5

##数据质控
scRNA3 <- subset(scRNA3, subset = nFeature_RNA > minGene)
col.num <- length(levels(as.factor(scRNA3@meta.data$orig.ident)))
violin <-VlnPlot(scRNA3, group.by = "orig.ident",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ercc"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)

violin=VlnPlot(object = scRNA3, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.ercc","percent.mt"), 
               group.by = 'orig.ident',
               ncol = 4) 
ggsave("QC/vlnplot_group_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_group_after_qc.png", plot = violin, width = 12, height = 6)

plot1 <- FeatureScatter(scRNA3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA3, feature1 = "nCount_RNA", feature2 = "percent.ercc")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("QC/pearplot_after_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot_after_qc.png", plot = pearplot, width = 12, height = 5)

