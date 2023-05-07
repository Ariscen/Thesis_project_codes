library(ggplot2)

eqtl_geno <- readRDS("eqtl_geno.rds")
compare_figure_celltype_ind <- readRDS("compare_figure_celltype_ind.rds")
simu_sce1 <- readRDS("simu_sce1.rds")
example_sce1 <- readRDS("example_sce1.rds")

plot_feature <- function(compare_figure_celltype_ind,example_sce1,simu_sce1,feature){
  compare_figure_celltype_ind$feature <- c(logcounts(example_sce1)[feature,],
                                           logcounts(simu_sce1)[feature,])
  ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "feature")) + 
    geom_point(alpha = 0.5, size = 1) + scale_color_continuous(type = "viridis") +
    facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"LYZ")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"CD8A")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"IL7R")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"GZMA")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"TCL1A")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"FCER1A")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"GNLY")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"S100B")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"GZMK")
plot_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,"CD27")

plot_all<-function(compare_figure_celltype_ind,color_by){
  ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = color_by)) + 
    geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + scale_color_discrete(type = "viridis") +
    facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_all(compare_figure_celltype_ind,"celltype")
plot_all(compare_figure_celltype_ind,"ind")

plot_ind <- function(compare_figure_celltype_ind,example_sce1,simu_sce1,inds){
  compare_figure_celltype_ind <- compare_figure_celltype_ind[which(compare_figure_celltype_ind$ind%in%inds),]
  ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "celltype")) + 
    geom_point(alpha = 0.5, size = 1) + scale_color_discrete(type = "viridis") +
    facet_wrap(ind~Method,nrow = length(inds)) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_ind(compare_figure_celltype_ind,example_sce1,simu_sce1,inds = c("1_1","2_2"))

plot_ind_feature <- function(compare_figure_celltype_ind,example_sce1,simu_sce1,feature,inds){
  compare_figure_celltype_ind$feature <- c(logcounts(example_sce1)[feature,],
                                           logcounts(simu_sce1)[feature,])
  compare_figure_celltype_ind <- compare_figure_celltype_ind[which(compare_figure_celltype_ind$ind%in%inds),]
  ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "feature")) + 
    geom_point(alpha = 0.5, size = 1) + scale_color_continuous(type = "viridis") +
    facet_wrap(ind~Method,nrow = length(inds)) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_ind_feature(compare_figure_celltype_ind,example_sce1,simu_sce1,
             feature = "CD8A",inds = c("1_1","2_2"))

set.seed(123)
compare_figure_combined <- plot_reduceddim(ref_sce = example_sce1, 
                                           sce_list = list(simu_sce_ref1,simu_sce1), 
                                           name_vec = c("Reference", "scDesign3","scDesign-pop"),
                                           assay_use = "logcounts", 
                                           if_plot = FALSE, 
                                           color_by = "celltype", 
                                           n_pc = 20)

compare_figure_combined$ind <- rep(colData(example_sce1)[,"ind"],3)
saveRDS(compare_figure_combined,"compare_figure_combined.rds")

plot_all(compare_figure_combined,"celltype")
plot_all(compare_figure_combined,"ind")

plot_ind(compare_figure_combined,example_sce1,simu_sce1,inds = c("1_1","2_2"))

plot_ind_celltype <- function(compare_figure_celltype_ind,example_sce1,simu_sce1,inds,celltypes,plot.type){
  compare_figure_celltype_ind <- compare_figure_celltype_ind[which(compare_figure_celltype_ind$ind%in%inds &
                                                                     compare_figure_celltype_ind$celltype%in%celltypes),]
  if(plot.type == "pca"){
    ggplot(compare_figure_celltype_ind, aes_string(x = "PC1", y = "PC2",color = "ind")) + 
      geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + 
      facet_wrap(ind~Method,ncol = 6) + theme(aspect.ratio = 1, legend.position = "bottom") + 
      guides(color = guide_legend(override.aes = list(size = 2,alpha = 1)))
  }
  else if(plot.type == "umap"){
    ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "ind")) + 
      geom_point(alpha = 0.5, size = 1) +
      facet_wrap(ind~Method,ncol = 6) + theme(aspect.ratio = 1, legend.position = "bottom") 
  }
}

plot_ind_celltype(compare_figure_celltype_ind = compare_figure_combined,
                  example_sce1 = example_sce1,
                  simu_sce1 = simu_sce1,
                  inds = unique(compare_figure_celltype_ind$ind)[1:9],
                  celltypes = unique(compare_figure_celltype_ind$celltype)[5],
                  plot = "umap")

### 
sce_combined <- SingleCellExperiment(
  assay = list(counts = cbind(counts(example_sce1), counts(simu_sce1)),
                              logcounts = cbind(logcounts(example_sce1), logcounts(simu_sce1)))
)

colnames(sce_combined) <- c(paste0(colnames(example_sce1),"-r"),paste0(colnames(simu_sce1),"-p"))

compare_figure_combined <- data.frame(compare_figure_combined)
rownames(compare_figure_combined) <- c(paste0(colnames(example_sce1),"-r"),
                                       paste0(colnames(example_sce1),"-o"),
                                       paste0(colnames(example_sce1),"-p"))

methods <- unique(compare_figure_combined$Method)
methods
colData(sce_combined)$Method <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,3)]),1]
colData(sce_combined)$celltype <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,3)]),2]
colData(sce_combined)$ind <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,3)]),25]

reducedDim(sce_combined, "PCA") <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,3)]),3:22]
reducedDim(sce_combined, "UMAP") <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,3)]),23:24]

sce_combined
saveRDS(sce_combined,"sce_combined.rds")
sce_combined <- readRDS("sce_combined.rds")

sce_pca <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined, 
                                     k = 50, n_dim = 2, cell_min = 4,
                                     res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "PCA")
sce_umap <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined, 
                                      k = 50, n_dim = 2, cell_min = 4,
                                      res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "UMAP")
PCA_lisi <- mean(colData(sce_pca)$weighted_isi)
PCA_entropy <- mean(colData(sce_pca)$entropy)
PCA_cms <- ks.test(colData(sce_pca)$cms.cms,"punif", 0, 1)$statistic

UMAP_lisi <- mean(colData(sce_umap)$weighted_isi)
UMAP_entropy <- mean(colData(sce_umap)$entropy)
UMAP_cms <- ks.test(colData(sce_umap)$cms.cms,"punif", 0, 1)$statistic

res <- c(PCA_lisi, PCA_entropy, PCA_cms, UMAP_lisi, UMAP_entropy, UMAP_cms)
names(res) <- c("PCA_lisi", "PCA_entropy", "PCA_cms", "UMAP_lisi", "UMAP_entropy", "UMAP_cms")
res

saveRDS(res,"res.rds")
saveRDS(sce_pca,"sce_pca.rds")
saveRDS(sce_umap,"sce_umap.rds")

inds <- sort(unique(colData(sce_combined)$ind))
inds
saveRDS(inds,"inds.rds")
ind_res <- pbapply::pblapply(inds,function(ind){
  sce_combined_sel <- sce_combined[,which(colData(sce_combined)$ind==ind)]
  sce_pca <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined_sel, 
                                       k = 50, n_dim = 2, cell_min = 4,
                                       res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "PCA")
  sce_umap <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined_sel, 
                                        k = 50, n_dim = 2, cell_min = 4,
                                        res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "UMAP")
  PCA_lisi <- mean(colData(sce_pca)$weighted_isi)
  PCA_entropy <- mean(colData(sce_pca)$entropy)
  PCA_cms <- ks.test(colData(sce_pca)$cms.cms,"punif", 0, 1)$statistic
  
  UMAP_lisi <- mean(colData(sce_umap)$weighted_isi)
  UMAP_entropy <- mean(colData(sce_umap)$entropy)
  UMAP_cms <- ks.test(colData(sce_umap)$cms.cms,"punif", 0, 1)$statistic
  
  res <- c(PCA_lisi, PCA_entropy, PCA_cms, UMAP_lisi, UMAP_entropy, UMAP_cms)
  names(res) <- c("PCA_lisi", "PCA_entropy", "PCA_cms", "UMAP_lisi", "UMAP_entropy", "UMAP_cms")
  return(res)
})

saveRDS(ind_res,"ind_res.rds")
