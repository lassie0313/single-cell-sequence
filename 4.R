pc.num = 1:30
sce_bone <- NormalizeData(sce_bone,verbose = FALSE)
sce_bone <- FindVariableFeatures(sce_bone,selection.method = "vst", nfeatures = 600,verbose = FALSE)
sce_bone <- ScaleData(sce_bone,verbose = FALSE)
sce_bone <- RunPCA(sce_bone,npcs = 30,verbose = FALSE)
sce_bone <- RunUMAP(sce_bone,reduction = "pca",dims = pc.num)
sce_bone = RunTSNE(sce_bone, dims = pc.num)
p1 <- DimPlot(sce_bone,reduction = "umap",group.by = "orig.ident")
plot(p1)
p2 <- DimPlot(sce_bone,reduction = "tsne",group.by = "orig.ident")
plot(p2)

scRNAlist <- SplitObject(sce_bone,split.by = "orig.ident")
for(i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]],verbose = FALSE)
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = "vst",nfeatures = 600,
                                         verbose = FALSE)
}
reference.list <- scRNAlist[c("E11","E13","E15","E18")]
scRNA.anchors <- FindIntegrationAnchors(object.list = reference.list,dims = 1:30)
scrna.integrated <- IntegrateData(anchorset = scRNA.anchors,dims = 1:30)
#切换到整合后的assay
DefaultAssay(scrna3) <- "integrated"
