
rm(list = ls())
options(stringsAsFactors = F)

# 加载原始表达矩阵
load(file = "../Analysis/data/Step01-airwayData.Rdata")

# 读取3个软件的差异分析结果
load(file = "../Analysis/deg_analysis/Step03-limma_voom_nrDEG.Rdata")
load(file = "../Analysis/deg_analysis/Step03-DESeq2_nrDEG.Rdata")
load(file = "../Analysis/deg_analysis/Step03-edgeR_nrDEG.Rdata")
ls()

# 根据需要修改DEG的值
fc_cutoff <- 1
fdr <- 50
data <- DEG
colnames(data)
data$regulated <- "normal"

loc_up <- intersect(which(data$avg_log2FC>log2(fc_cutoff)),which(data$FDR>fdr))
loc_down <- intersect(which(data$avg_log2FC< (-log2(fc_cutoff))),which(data$FDR>fdr))

data$regulated[loc_up] <- "up"
data$regulated[loc_down] <- "down"

table(data$regulated)
data$FDR=-log10(data$p_val_adj)


data2=data[which(data$FDR<150),]
f=which(!(data$p_val_adj==0))
data3
# 绘制火山图
library(ggplot2)
library(ggrepel)
colnames(data)
label=c("Cxcr4", "Osm","Tgfb1","Vegfa","Adgre5","Cd83","Sirpa","Nrg1")
p <- ggplot(data=data2, aes(x=avg_log2FC, y=-log10(p_val_adj),color=regulated)) + 
     geom_point(alpha=0.5, size=1.8) + theme_set(theme_set(theme_bw(base_size=20))) + 
     xlab("log2FC") + ylab("-log10(FDR)") +scale_colour_manual(values = c('blue','black','red'))
p+geom_text_repel(data=data2,aes(x=avg_log2FC, y=-log10(p_val_adj),label=label),size=2)
p







