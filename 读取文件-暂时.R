require(openxlsx)
file_names <- dir("F:/单细胞测序/老板/DEMO1/新建文件夹/Marrow")
pattern = "*.csv", recursive = F, full.names = T)
# 创建存储数据的数据框，直接将第一个文件的数据赋值给它
df <- read.csv(file_names[1])
# 从第二个文件开始合并
for (i in 2:length(file_names)) {
  df <- rbind(df, read.csv(file_names[i]))}
a=read.csv("新建文件夹/Marrow/A1-D041912-3_8_M-1-1.csv")
