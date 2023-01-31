rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
library(monocle)
library(Seurat)
options(stringsAsFactors=F)

#其他目标细胞类型并保存（scRNA为之前构建的Seurat对象）
#tcell <- subset(scRNA, idents = 
                  #c("Naive CD4 T","Memory CD4 T","CD8 T"))
#saveRDS(tcell,"tcell.rds")

scRNA <- readRDS("scRNA.rds")
str(scRNA)
## 第1步：构建monocle对象
##提取原始的表达矩阵：UMI count
expr_matrix <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
head(expr_matrix[,1:4])
dim(expr_matrix)

##提取表型信息--细胞信息
p_data <- scRNA@meta.data
head(scRNA@active.ident)
##整合每个细胞的细胞鉴定信息到p_data里面
p_data$celltype <- scRNA@active.ident
head(p_data)
dim(p_data)

##提取基因信息
f_data <- data.frame(gene_short_name = row.names(scRNA), row.names = row.names(scRNA))
head(f_data)

#expr_matrix的行数与featureData的行数相同（gene number）
#expr_matrix的列数与phenoData的行数相同(cell number)
ncol(expr_matrix)
nrow(f_data)
nrow(expr_matrix)
nrow(p_data)
##构建monocle2对象
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,      #表达矩阵信息（10x为UMI count值）
                      phenoData = pd,   #表型信息，主要为每个细胞的信息，建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息
                      featureData = fd, #基因信息，如生物类型、gc含量等
                      lowerDetectionLimit = 0.5, #表达量的下限
                      expressionFamily = negbinomial.size()) #数据的分布（不同的数据类型有不同的分布），适用于不同的统计检验方法

#expressionFamily参数说明：
#1.negbinomial.size()和negbinomial()：输入的表达矩阵为UMI,一般适用于10x的数据；negbinomial()的结果更准确，但是计算更耗时；一般建议采用negbinomial.size()。
#2.tobit():适用于输入的表达矩阵为FPKM或者TPM, 构建monocle2的class时会自动进行log化计算
#3.gaussianff():输入为log化后的FPKM或者TPM
#目前在单细胞数据中，FPKM已不多用，smart-seq2平台数据一般采用TPM


## 第2步：估计文库大小及分散度 (类似于seurat的归一化处理）
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## 第3步：细胞质控，因为我们是基于Seurat的结果进行后续分析，Seurat已经完成了细胞质控和过滤，因此省略


## 第4步：选择用来定义轨迹的基因及可视化
#选择有意义的基因集，来构建细胞轨迹（类似与细胞聚类中的高变基因选择进行细胞聚类）
#基因集的选择是轨迹构建结果的最决定因素

#为了加速计算，我们首先过滤掉表达的细胞数过少的基因

##计算每个基因表达的细胞数目
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

##挑选表达5个及以上的细胞数的基因list
express_genes <- row.names(subset(fData(cds),num_cells_expressed>=5))
head(express_genes)
##对符合要求的基因做差异分析
diff <- differentialGeneTest(cds[express_genes,],fullModelFormulaStr="~celltype",cores=2)
head(diff)
##如果想分析其他的改~celltype
#diff_Estatus <- differentialGeneTest(cds[express_genes,],fullModelFormulaStr="~orig.ident",cores=2)
#head(diff_Estatus)

##差异表达基因作为轨迹构建的基因
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

deg_Estatus <- subset(diff_Estatus, qval < 0.01)
deg_Estatus <- deg_Estatus[order(deg_Estatus$qval,decreasing=F),]
head(deg_Estatus)
##差异基因的结果文件保存
write.table(deg,file="scRNA.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
write.table(deg_Estatus,file="scRNA.monocle.DEG_Estatus.xls",col.names=T,row.names=F,sep="\t",quote=F)

####注：第52-59行可能运行比较耗时，若为了节省时间，可以直接载入已经生成的差异基因列表，代码如下：
#deg <- read.table("../Data/scRNA.monocle.DEG.xls",header=T,sep="\t",check.names=F)

#基因集选择的几种常见策略：
#1）基于Seurat 鉴定的marker genes
#2）基于monocle2差异分析鉴定的差异基因，这次的属于这个
#3）基于高变基因
#4）其他自定义的基因集
#---基因的选择比较灵活，但是要有生物学或者数据的依据

## 轨迹构建基因可视化
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
pdf("scRNA.ordergenes.pdf")
 plot_ordering_genes(cds)
dev.off()
ggsave("scRNA.ordergenes.png",plot_ordering_genes(cds))

ordergene_Estatus <- rownames(deg_Estatus)
cds <- setOrderingFilter(cds, ordergene)
pdf("scRNA.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
ggsave("scRNA.ordergenes.png",plot_ordering_genes(cds))

## 第5步：数据降维
cds <- reduceDimension(cds,max_components = 2, method = 'DDRTree')

## 第6步：将细胞安装拟时进行排序，构建轨迹
cds <- orderCells(cds)

## 第7步：轨迹可视化
##以pseudotime值上色  Pseudotime是monocle2基于细胞基因表达信息计算的概率，表示时间的先后。

pdf("scRNA.monocle.pseudotime.pdf")
 plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
dev.off()
ggsave("scRNA.monocle.pseudotime.png",plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE))

##以细胞类型上色
pdf("scRNA.monocle.cluster.pdf")
 plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE)
dev.off()
ggsave("scRNA.monocle.cluster.png",plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE))

##以细胞状态上色
pdf("scRNA.monocle.state.pdf")
 plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()
ggsave("scRNA.monocle.state.png",plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE))

##以细胞状态上色（拆分）
pdf("scRNA.monocle.state.faceted.pdf")
 plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
dev.off()
ggsave("scRNA.monocle.state.faceted.png",plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1))

## 第8步 :关键驱动基因的表达变化图
#！！关键基因的确定需要结合项目背景，应为在轨迹分化中起关键作用的基因
#选择前4个 top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]

pdf("scRNA.keygene.state.pdf")
 plot_genes_in_pseudotime(cds_subset, color_by = "State")
dev.off()
ggsave("scRNA.keygene.pseu.png",plot_genes_in_pseudotime(cds_subset, color_by = "State"))

pdf("scRNA.keygene.celltype.pdf")
 plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
dev.off()
ggsave("scRNA.keygene.celltype.png",plot_genes_in_pseudotime(cds_subset, color_by = "celltype"))

pdf("scRNA.keygene.pseudotime.pdf")
 plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
dev.off()
ggsave("scRNA.keygene.pseudotime.png",plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime"))

## 第9步：基于表达趋势相似度的基因聚类，选择前30个差异基因
topgene <- head(ordergene,30)
pdf("scRNA.pheatmap.pdf")
 plot_pseudotime_heatmap(cds[topgene,],num_clusters = 3,cores = 1,show_rownames = T)
dev.off()
ggsave("scRNA.pheatmap.png")

## 第10步：分支分析
#注：第10步与前面9步相对独立，进行此分析，需明确细胞确实存在这分支分化，即同一细胞类型转为2类或以上
##BRAM函数输入的是经过order后的轨迹结果对象
BEAM_res <- BEAM(cds, branch_point = 2, cores = 1) #表示有三个分支，根据之前的轨迹图可知
BEAM_res <- BEAM_res[order(BEAM_res$qval,decreasing=F),] #基于qval值进行排序（从小到大）
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")] #选择特定的3列
head(BEAM_res)

##分支相关基因热图展示 最好的是中间灰色为原始的细胞，有高表达的基因 左右两边cell line都有高表达的基因
pdf("scRNA.branched_heatmap.pdf")
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                          qval < 1e-4)),], #表示选择此阈值下的差异基因
                                          branch_point = 1, #表示有3个分支
                                          num_clusters = 4, #表示基因聚为4类（可以调整）
                                          cores = 2,
                                          use_gene_short_name = T, #表示是否采用gene name
                                          show_rownames = T) #表示是否展示基因名称（或者geneid）
dev.off()
ggsave("scRNA.branched_heatmap.png",plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],branch_point = 3,num_clusters = 4,cores = 2,use_gene_short_name = T,show_rownames = T))

##关键基因的表达分支趋势图
#选择前4个差异基因
branch_gene <- rownames(BEAM_res[1:4,])
#以celltype进行上色和可视化，其余两种上色方案类似，不再详述
pdf("scRNA.branched_genes.pdf")
 plot_genes_branched_pseudotime(cds[branch_gene,],
                        branch_point = 1,
                        color_by = "celltype",
                        ncol = 1)

dev.off()
ggsave("scRNA.branched_genes.png",plot_genes_branched_pseudotime(cds[branch_gene,],branch_point = 1,color_by = "celltype",ncol = 1))

##保存monocle2对象
saveRDS(cds,"scRNA.Tcell.monocle2.rds")
