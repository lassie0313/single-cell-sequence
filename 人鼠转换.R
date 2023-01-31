m_gene=read.table("m_gene.txt")
library(biomaRt)

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = m_gene , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
genes <- convertMouseGeneList(m_gene)

gene=read.table("gene.txt")
counts_final=counts[match(gene$MGI.symbol,counts$Gene),]
a=as.data.frame(a)
b=gene[match(m_)]
