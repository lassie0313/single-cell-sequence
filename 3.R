install.packages('devtools')
.libPaths()
packageVersion("Seurat")
library(Seurat)

Sys.setenv(R_MAX_NUM_DLLS=999)
samples=list.files("wt/")
samples
dir <- file.path('./wt/',samples)
names(dir) <- samples
scRNAlist <- list()
for(i in 1:length(dir)){
  print(i)
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.genes=1000)
}
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], 
                                    scRNAlist[[4]]))
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数

library(Seurat)
# 步骤 ScaleData 的耗时取决于电脑系统配置（保守估计大于一分钟）
scRNA2 <- ScaleData(object = scRNA2, 
                     vars.to.regress = c('nCount_RNA'), 
                     model.use = 'linear', 
                     use.umi = FALSE)
scRNA2 <- FindVariableFeatures(object = scRNA2, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, 
                                x.high.cutoff = 4, 
                                y.cutoff = 0.5)
length(VariableFeatures(scRNA2)) 
scRNA2 <- RunPCA(object = scRNA2, pc.genes = VariableFeatures(scRNA2))
# 下面只是展现不同降维算法而已，并不要求都使用一次。
scRNA2 <- RunICA(sce.big )
scRNA2 <- RunTSNE(sce.big )
#sce.big <- RunUMAP(sce.big,dims = 1:10)
#VizPCA( sce.big, pcs.use = 1:2)
DimPlot(object = sce.big, reduction = "pca") 
DimPlot(object = sce.big, reduction = "ica")
DimPlot(object = sce.big, reduction = "tsne")