##都在linux上做
library(dplyr)
scrna=readRDS("scRNA3_recombine.rds")
write.table(data.frame(cell=rownames(scrna@meta.data),celltype=scrna@meta.data$celltype),'meta.txt',sep='\t',quote=F,row.names=F)
write.table(dplyr::bind_cols(Gene=rownames(scrna@assays$RNA@data),as.data.frame(scrna@assays$RNA@data)),'counts.txt',sep='\t',quote=F,row.names=F)

###cellphone db只能人 用biomart包将鼠变成人
library(data.table)
counts=fread("counts.txt")
m_gene=counts$Gene
library(biomaRt)
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = m_gene , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
write.table(genesV2,file="gene.txt") 
gene=read.table("gene.txt")
counts_final=counts[match(gene$MGI.symbol,counts$Gene),]
counts_final[1:5,1:5]
#     Gene S1_AAACCCATCCTCGCAT-1 S1_AAAGTCCGTATGTCCA-1 S1_AAGCGTTTCTCAGGCG-1
#1: mt-Nd2              4.664316              3.708814               0.00000
#2: mt-Nd3              2.843056              2.287029               0.00000
#3: mt-Nd5              0.000000              1.690616               0.00000
#4: mt-Co3              5.899304              4.674208               5.08459
#5: mt-Nd4              4.737750              4.273399               2.83623
#S1_AAGGAATCATGCAGCC-1
#1:              4.135025
#2:              0.000000
#3:              3.180517
#4:              4.937025
#5:              4.687759

gene[1:5,1:2]
#MGI.symbol HGNC.symbol
#1     mt-Nd2      MT-ND2
#2     mt-Nd3      MT-ND3
#3     mt-Nd5      MT-ND5
#4     mt-Co3      MT-CO3
#5     mt-Nd4      MT-ND4

identical(counts_final$Gene,gene$MGI.symbol)
counts_final$Gene=gene$HGNC.symbol
counts_final[1:5,1:5]
#Gene S1_AAACCCATCCTCGCAT-1 S1_AAAGTCCGTATGTCCA-1 S1_AAGCGTTTCTCAGGCG-1
#1: MT-ND2              4.664316              3.708814               0.00000
#2: MT-ND3              2.843056              2.287029               0.00000
#3: MT-ND5              0.000000              1.690616               0.00000
#4: MT-CO3              5.899304              4.674208               5.08459
#5: MT-ND4              4.737750              4.273399               2.83623
#S1_AAGGAATCATGCAGCC-1
#1:              4.135025
#2:              0.000000
#3:              3.180517
#4:              4.937025
#5:              4.687759
write.table(counts_final,file="counts_final.txt",sep='\t',quote=F,row.names=F)

#import to linux 退出R
cellphonedb method statistical_analysis meta.txt counts_final.txt --counts-data=gene_name
cellphonedb plot dot_plot --means-path /home/data/gma04/scrna/marrow/monocle/out/means.txt --pvalues-path /home/data/gma04/scrna/marrow/monocle/out/pvalues.txt --output-name Dot_Plot——msc.pdf --rows row.txt 
cellphonedb plot heatmap_plot  meta.txt  --pvalues-path /home/data/gma04/scrna/marrow/monocle/out/pvalues.txt  --count-name Heatmap_Count.pdf --log-name Heatmap_Log_count.pdf --count-network-name count_network.txt

install.packages("psych")
install.packages("qgraph")
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
mynet<- read.delim("count_network.txt", check.names = FALSE)
mynet %>% filter(count>0) -> mynet
net   <- graph_from_data_frame(mynet)
print(net)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
dpi           <- 300
karate_groups <- cluster_optimal(net)
coords        <- layout_in_circle(net, order = order(membership(karate_groups)))  # 设置网络布局
E(net)$width  <- E(net)$count/10
net2          <- net
for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
}


png('Global_Network.png',w=12*dpi,h=8*dpi,units = "px",res = dpi,type='cairo')
p1 <- plot(net, edge.arrow.size=.1,
           edge.curved=0.2,
           vertex.color=allcolour,
           vertex.frame.color="#555555",
           vertex.label.color="black",
           layout = coords,
           vertex.label.cex=.7)
dev.off()

png('Global_Shell_Network.png',w=12*dpi,h=8*dpi,units = "px",res = dpi,type='cairo')
length(unique(mynet$SOURCE)) 
par(mfrow=c(2,5), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       edge.label = E(net1)$count, # 绘制边的权重
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1
  ) 
  
}

dev.off()

