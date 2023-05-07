library(SingleCellExperiment)
setwd("C:/Future_Dpan/sc-eQTL/data_backup/")
compare_figure_combined <- readRDS("compare_figure_combined.rds")
example_sce1 <- readRDS("example_sce1.rds")
simu_sce_ref1 <- readRDS("simu_sce_ref1.rds")

sce_combined_ref <- SingleCellExperiment(
  assay = list(counts = cbind(counts(example_sce1), counts(simu_sce_ref1)),
               logcounts = cbind(logcounts(example_sce1), logcounts(simu_sce_ref1)))
)

colnames(sce_combined_ref) <- c(paste0(colnames(example_sce1),"-r"),paste0(colnames(simu_sce_ref1),"-o"))

compare_figure_combined <- data.frame(compare_figure_combined)
rownames(compare_figure_combined) <- c(paste0(colnames(example_sce1),"-r"),
                                       paste0(colnames(example_sce1),"-o"),
                                       paste0(colnames(example_sce1),"-p"))

methods <- unique(compare_figure_combined$Method)
methods
colData(sce_combined_ref)$Method <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,2)]),1]
colData(sce_combined_ref)$celltype <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,2)]),2]
colData(sce_combined_ref)$ind <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,2)]),25]

reducedDim(sce_combined_ref, "PCA") <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,2)]),3:22]
reducedDim(sce_combined_ref, "UMAP") <- compare_figure_combined[which(compare_figure_combined$Method%in%methods[c(1,2)]),23:24]

sce_combined_ref
saveRDS(sce_combined_ref,"sce_combined_ref.rds")

sce_combined_ref <- readRDS("sce_combined_ref.rds")

sce_pca_ref <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined_ref, 
                                     k = 50, n_dim = 2, cell_min = 4,
                                     res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "PCA")
sce_umap_ref <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined_ref, 
                                      k = 50, n_dim = 2, cell_min = 4,
                                      res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "UMAP")
PCA_lisi_ref <- mean(colData(sce_pca_ref)$weighted_isi)
PCA_entropy_ref <- mean(colData(sce_pca_ref)$entropy)
PCA_cms_ref <- ks.test(colData(sce_pca_ref)$cms.cms,"punif", 0, 1)$statistic

UMAP_lisi_ref <- mean(colData(sce_umap_ref)$weighted_isi)
UMAP_entropy_ref <- mean(colData(sce_umap_ref)$entropy)
UMAP_cms_ref <- ks.test(colData(sce_umap_ref)$cms.cms,"punif", 0, 1)$statistic

res_ref <- c(PCA_lisi_ref, PCA_entropy_ref, PCA_cms_ref, UMAP_lisi_ref, UMAP_entropy_ref, UMAP_cms_ref)
names(res_ref) <- c("PCA_lisi", "PCA_entropy", "PCA_cms", "UMAP_lisi", "UMAP_entropy", "UMAP_cms")
res_ref

saveRDS(res_ref,"res_ref.rds")
saveRDS(sce_pca_ref,"sce_pca_ref.rds")
saveRDS(sce_umap_ref,"sce_umap_ref.rds")

inds <- sort(unique(colData(sce_combined_ref)$ind))
inds
ind_res_ref <- pbapply::pblapply(inds,function(ind){
  sce_combined_ref_sel <- sce_combined_ref[,which(colData(sce_combined_ref)$ind==ind)]
  sce_pca <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined_ref_sel, 
                                       k = 50, n_dim = 2, cell_min = 4,
                                       res_name = c("weighted_isi", "entropy", "cms"), group = "Method", dim_red = "PCA")
  sce_umap <- CellMixS::evalIntegration(metrics = c("isi", "entropy", "cms"), sce_combined_ref_sel, 
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

saveRDS(ind_res_ref,"ind_res_ref.rds")
