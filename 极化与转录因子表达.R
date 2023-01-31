#final
scRNA_macro_two=readRDS("scRNA_macro_two.rds")
a=c("Cebpb", "Klf3","Atf3","Atf4", "Fosl2","Klf4")
m=c("Il10")
p3=FeaturePlot(scRNA_macro_two, features =a, split.by = "position", max.cutoff = 3, 
               cols = c("grey", "red"))
plots <- VlnPlot(scRNA_macro_two, features = a, split.by = "position", group.by = "position", 
                 pt.size = 0, combine = FALSE)
plots=wrap_plots(plots = plots, ncol = 1)
#Mrc1
#b=c("Arg1","Cd86","Tnf","Nos2")
b=c("Nos2")

p3=FeaturePlot(scRNA_macro_two, features =b, split.by = "position", max.cutoff = 3, 
               cols = c("grey", "red"))
plots <- VlnPlot(scRNA_macro_two, features = b, split.by = "position", group.by = "position", 
                 pt.size = 0, combine = FALSE)
plots=wrap_plots(plots = plots, ncol = 1)


b=c("Nos2")
p4=FeaturePlot(scRNA_macro_two, features =b, split.by = "position", max.cutoff = 3, 
               cols = c("grey", "red"))
plots <- VlnPlot(scRNA_macro_two, features = b, split.by = "position", group.by = "position", 
                 pt.size = 0, combine = FALSE)
plots=wrap_plots(plots = plots, ncol = 1)
saveRDS(scRNA_macro_two,file="scRNA_macro_two.rds")


#M1极化 
#cd86 tnfa(Tnf) inos(Nos2) 

#M2极化
#Arg1 cd206(mmr mrc1)