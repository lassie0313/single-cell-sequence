
setwd("F:/单细胞测序/老板/DEMO1/longbone/Marrow")
files = list.files()
write.csv(files,file="filename.csv")
time1 = Sys.time()
for(i in 1:length(files)){
  if(i == 1) {
    temp = read.csv(files[i],header = F)
    name <- substr(files[i], start = 1, stop = 14)
    names(temp)[2] <- name
  } else {
    dt <- read.csv(files[i],header = F)
    name <- substr(files[i], start = 1, stop = 14)
    names(dt)[2] <- name
    temp = merge(temp, dt, by='V1')  }}

time2 = Sys.time()
time2 - time1
namere=read.csv("../filename.csv")
a=namere$x
colnames(temp)=a
rownames(temp)=temp$Gene
temp1=temp[,-1]
temp1[1:5,1:5]
library(Seurat)
scRNA3=readRDS("scRNA3.rds")
#scRNA3 = CreateSeuratObject(temp1)
scRNA3[["percent.ercc"]] <- PercentageFeatureSet(scRNA3, pattern = "^Ercc")
