library(ggplot2)
library(SingleCellExperiment)
setwd("C:/Future_Dpan/sc-eQTL/data_backup/")

eqtl_geno <- readRDS("eqtl_geno.rds")
compare_figure_combined <- readRDS("compare_figure_combined.rds")
simu_sce1 <- readRDS("simu_sce1.rds")
simu_sce_ref1 <- readRDS("simu_sce_ref1.rds")
example_sce1 <- readRDS("example_sce1.rds")

plot_feature <- function(compare_figure,example_sce1,simu_sce_ref1,simu_sce1,feature){
  compare_figure$feature <- c(logcounts(example_sce1)[feature,],
                              logcounts(simu_sce_ref1)[feature,],
                              logcounts(simu_sce1)[feature,])
  print(cor(logcounts(example_sce1)[feature,],logcounts(simu_sce_ref1)[feature,],method="pearson"))
  print(cor(logcounts(example_sce1)[feature,],logcounts(simu_sce1)[feature,],method="pearson"))
  #
  compare_figure <- compare_figure[which(compare_figure$Method!="scDesign3"),]
  #
  ggplot(compare_figure, aes_string(x = "UMAP1", y = "UMAP2",color = "feature")) + 
    geom_point(alpha = 0.5, size = 1) + scale_color_continuous(type = "viridis") +
    facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"LYZ")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"CD8A")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"IL7R")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"GZMA")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"TCL1A")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"FCER1A")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"GNLY")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"S100B")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"GZMK")
plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"CD27")

print(cor(as.numeric(as.vector(logcounts(example_sce1))),as.numeric(as.vector(logcounts(simu_sce_ref1))),method="pearson"))
print(cor(as.numeric(as.vector(logcounts(example_sce1))),as.numeric(as.vector(logcounts(simu_sce1))),method="pearson"))
