
# 清空当前对象
rm(list = ls())
options(stringsAsFactors = F)

# 读取基因表达矩阵
lname <- load("C:/Users/lassie/Desktop/生信/转录组测序/code-down/code-down/data/Step01-airwayData.Rdata")

lname

exprSet <- filter_count #pca 是用的标准化后的值,这边用的是原始的 count 值
# 检查表达谱
dim(exprSet)
exprSet[1:6,1:6]
table(group_list)

# 加载包
library(limma)
library(edgeR)

## 第一步，创建设计矩阵和对比：假设数据符合正态分布，构建线性模型
# 0代表x线性模型的截距为0
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(exprSet)
design

# 设置需要进行对比的分组，需要修改
comp <- 'trt-untrt' #哪个写在前面哪个就是实验组
cont.matrix <- makeContrasts(contrasts=c(comp),levels = design)



## 第二步，进行差异表达分析
# 将表达矩阵转换为edgeR的DGEList对象
dge <- DGEList(counts=exprSet)

# 进行标准化
dge <- calcNormFactors(dge)   

#Use voom() [15] to convert the read counts to log2-cpm, with associated weights, ready for linear modelling:
v <- voom(dge,design,plot=TRUE, normalize="quantile") 
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)


## 第三步，提取过滤差异分析结果
tmp <- topTable(fit2, coef=comp, n=Inf,adjust.method="BH")
DEG_limma_voom <- na.omit(tmp)
head(DEG_limma_voom)

# 筛选上下调，设定阈值 看你想研究的基因在不在里面，如果不在可以放宽条件 log2（1.5）=0.58
fc_cutoff <- 2
fdr <- 0.05

DEG_limma_voom$regulated <- "normal"

loc_up <- intersect(which(DEG_limma_voom$logFC>log2(fc_cutoff)),which(DEG_limma_voom$adj.P.Val<fdr))
loc_down <- intersect(which(DEG_limma_voom$logFC< (-log2(fc_cutoff))),which(DEG_limma_voom$adj.P.Val<fdr))

DEG_limma_voom$regulated[loc_up] <- "up"
DEG_limma_voom$regulated[loc_down] <- "down"
  
table(DEG_limma_voom$regulated)

# 添加一列gene symbol
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_limma_voom), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
head(id2symbol)

symbol <- rep("NA",time=nrow(DEG_limma_voom))
symbol[match(id2symbol[,1],rownames(DEG_limma_voom))] <- id2symbol[,2]
DEG_limma_voom <- cbind(rownames(DEG_limma_voom),symbol,DEG_limma_voom)
colnames(DEG_limma_voom)[1] <- "GeneID"


# 保存
write.table(DEG_limma_voom,"../code-down/data/DEG_limma_voom_all-1.xls",row.names = F,sep="\t",quote = F)

## 取表达差异倍数和p值,矫正后的pvalue
DEG_limma_voom <- DEG_limma_voom[,c(1,2,3,6,7,9)]
save(DEG_limma_voom, file = "../code-down/data/Step03-limma_voom_nrDEG.Rdata")


## 检查是否上下调设置错了
# 挑选一个差异表达基因
head(DEG_limma_voom)

exp <- c(t(express_cpm[match("ENSG00000178695",rownames(express_cpm)),]))
test <- data.frame(value=exp,group=group_list)
library(ggplot2)
ggplot(data=test,aes(x=group,y=value,fill=group)) + geom_boxplot()




