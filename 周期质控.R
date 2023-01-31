#细胞周期评分
mouse.cc.genes=readRDS("mouse_cell_cycle_genes.rds")
m.s.genes <- mouse.cc.genes
m.g2m.genes <- mouse.cc.genes


s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce_E11=CellCycleScoring(object = sce_E11, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce_E11, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0.1)
ggsave(filename="Vlnplot4_cycle.pdf",plot=p4)

sce_E11@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
ggsave(filename="cycle_details.pdf" )
