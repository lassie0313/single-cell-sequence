

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GENIE3", "AUCell", "RcisTarget"), version = "3.12")
install.packages('zoo')
BiocManager::install(c("mixtools", "rbokeh"))
BiocManager::install(c("NMF", "pheatmap", "Rtsne", "R2HTML"))
BiocManager::install(c("doMC", "doRNG"))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
BiocManager::install(c("SingleCellExperiment"))
devtools::install_github("aertslab/SCENIC", ref="v1.1.0")
packageVersion("SCENIC")
library(SCENIC)
library(Seurat)
library(RcisTarget)
library(GENIE3)
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather","https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed

,"https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather"
feather_database_url='https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather'
feather_database="${feather_database_url##*/}"
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL))
  descrURL <- gsub(".feather$", ".descr", featherURL)
  if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
}

awk -v feather_database=${feather_database} '$2 == feather_database' sha256sum.txt | sha256sum -c -


cellInfo <- data.frame(scRNA_macro_two@meta.data)
view(scRNA_macro_two@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="cellInfo.Rds")
##准备表达矩阵
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA_macro_two),1000)
exprMat1 <- scRNA_macro_two[,subcell]
DimPlot(exprMat1, reduction = "tsne", label = T,group.by = 'position')
table(exprMat1@meta.data$position)
saveRDS(exprMat1, "exprMat1.rds")
exprMat1=readRDS("exprMat1.rds")

install.packages("SCENIC",repos='https://mran.microsoft.com/snapshot/2019-02-01/')
library(SCENIC)
rm(list=ls())
dir.create("SCENIC")
dir.create("SCENIC/int")
setwd("~/SCENIC")

##准备scenic输入文件
#exprMat <- GetAssayData(exprMat1, assay = 'RNA', slot = 'data') %>% as.matrix()
exprMat <- as.matrix(exprMat1@assays$RNA@counts)

mydbDIR <- '/home/data/gma04/database/cisTarget/db'
dir(mydbDIR)
mydbs <- c('mm9-500bp-upstream-7species.mc9nr.feather','mm9-tss-centered-10kb-7species.mc9nr.feather')
names(mydbs) <- c("500bp", "10kb")
#小鼠org="mgi"
scenicOptions <- initializeScenic(org = "mgi", nCores = 16, dbDir = mydbDIR, datasetTitle = "macro")
saveRDS(scenicOptions, "int/scenicOptions.rds")
#scenicOptions = readRDS("int/scenicOptions.rds")
#scenicOptions@settings$nCores <- 10

##如果需要高性能计算服务，请保存以下文件联系Kinesin


##基因过滤
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 0.015 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]

##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
##计算TF-targets相关性
runGenie3(exprMat_filtered_log, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)

save(exprMat_filtered_log,scenicOptions,file = "input_GENIE3_data.Rdata")
##推断转录调控网络（regulon）
#此步运行时间长且极耗内存，慎用多线程！！！
scenicOptions@settings$nCores <- 2
library(BiocParallel)
runSCENIC_2_createRegulons(scenicOptions)

exprMat_all <- as.matrix(scRNA_macro_two@assays$RNA@data)
exprMat_all <- log2(exprMat_all+1)
scenicOptions@settings$nCores <- 1
library(doParallel)
scenicOptions <- initializeScenic(org="mgi", 
                                  dbDir=db , nCores=1) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)
#电脑内存只有16G，以上步骤相当耗时，可以睡前运行，醒来结果就出来了=。=!
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#网络活性的结果进行二维化
runSCENIC_4_aucell_binarize(scenicOptions)

scenicOptions=readRDS("scenicOptions.Rds")
#可视化
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))
## [1] "int/tSNE_AUC_05pcs_15perpl.Rds"

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_all, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Cebpb","Klf3","Atf3", "Atf4","Fosl2","Klf4")],], plots="Expression")

# Save AUC as PDF:
#linux
Cairo::CairoPDF("two/SCENIC/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()


#绘制密度图，显示稳定状态的细胞
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)
#同时显示多个调控元件
#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c("Cebpb","Klf3_extended","Klf4")
                  "Atf3", "Atf4_extended","Fosl2")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)

regulonNames <- list(red=c("Cebpb","Fosl2"),
                    green=c("Klf3_extended","Klf4"),
                    blue=c("Atf3","Atf4_extended"))

cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

#2.GRN:调控元件的靶点和motif
#查看regulon包含哪些基因
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Atf4_extended")]
 regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
#查看TF-target pair
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Cebpb" & highConfAnnot==TRUE]
#查看motif富集分析
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Cebpb"]
#3.查看cluster或者已知细胞类型的调控因子活性
scRNA_macro_two=readRDS("scRNA_macro_two.rds")
cellInfo <- data.frame(seuratCluster=Idents(scRNA_macro_two))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
topRegulators <- reshape2::melt(apply(regulonActivity_byCellType_Scaled, 2, function(x) cbind(sort(x[x>0], decreasing=TRUE))))[c("L1","Var1", "value")]; colnames(topRegulators) <- c("CellType","Regulon", "RelativeActivity")

#用二进制结果表示调控因子活性
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster),
                                              function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
p1=pheatmap::pheatmap(binaryActPerc_subset, color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA)

topRegulators <- reshape2::melt(apply(regulonActivity_byCellType_Binarized, 2, function(x) cbind(sort(x[x>minPerc], decreasing=TRUE))))[c("L1","Var1", "value")]; colnames(topRegulators) <- c("CellType","Regulon", "RelativeActivity")

#查看其它方法计算embeddings/trajectorie上的调控因子活性
dr_coords <- Embeddings(scRNA_macro_two, reduction="tsne")
tfs <- c("Fosl2")


par(mfrow=c(2,2))

load("cellsAUC.Rdata")#在liunx得到的
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=cellsAUC, plots = "AUC")
