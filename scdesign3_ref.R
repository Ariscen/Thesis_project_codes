library(scDesign3)
library(ggplot2)
setwd("C:/Future_Dpan/sc-eQTL/data_backup/")

example_sce1 <- readRDS("example_sce1.rds")

set.seed(123)
# example_simu_ref1 <- scdesign3(
#   sce = example_sce1,
#   assay_use = "counts",
#   celltype = "celltype",
#   pseudotime = NULL,
#   spatial = NULL,
#   other_covariates = NULL,
#   mu_formula = "celltype",
#   sigma_formula = "1",
#   family_use = "nb",
#   n_cores = 1,
#   usebam = FALSE,
#   corr_formula = "celltype",
#   copula = "gaussian",
#   DT = TRUE,
#   pseudo_obs = FALSE,
#   return_model = FALSE,
#   nonzerovar = FALSE,
#   parallelization = "pbmcmapply"
# )
colData(example_sce1)[,"group"] <- as.numeric(colData(example_sce1)[,"ind"])
Onek1k_data_ref1 <- construct_data(
  sce = example_sce1,
  assay_use = "counts",
  celltype = "celltype",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("ind"),
  corr_by = "group"
)

saveRDS(Onek1k_data_ref1,"Onek1k_data_ref1.rds")

Onek1k_marginal_ref1 <- fit_marginal(
  data = Onek1k_data_ref1,
  predictor = "gene",
  mu_formula = "celltype",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 1,
  usebam = FALSE,
  parallelization = "pbmcmapply"
)

saveRDS(Onek1k_marginal_ref1,"Onek1k_marginal_ref1.rds")

Onek1k_copula_ref1 <- fit_copula(
  sce = example_sce1,
  assay_use = "counts",
  marginal_list = Onek1k_marginal_ref1,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 2,
  new_covariate = NULL,
  input_data = Onek1k_data_ref1$dat
)

saveRDS(Onek1k_copula_ref1,"Onek1k_copula_ref1.rds")

Onek1k_para_ref1 <- extract_para(
  sce = example_sce1,
  marginal_list = Onek1k_marginal_ref1,
  n_cores = 1,
  family_use = "nb",
  new_covariate = NULL
)

saveRDS(Onek1k_para_ref1,"Onek1k_para_ref1.rds")

Onek1k_newcount_ref1 <- simu_new(
  sce = example_sce1,
  mean_mat = Onek1k_para_ref1$mean_mat,
  sigma_mat = Onek1k_para_ref1$sigma_mat,
  zero_mat = Onek1k_para_ref1$zero_mat,
  quantile_mat = NULL,
  copula_list = Onek1k_copula_ref1$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = Onek1k_data_ref1$dat,
  new_covariate = Onek1k_data_ref1$new_covariate,
  important_feature = Onek1k_copula_ref1$important_feature
)

saveRDS(Onek1k_newcount_ref1,"Onek1k_newcount_ref1.rds")

logcounts(example_sce1) <- log1p(counts(example_sce1))
simu_sce_ref1 <- example_sce1
counts(simu_sce_ref1) <- Onek1k_newcount_ref1
logcounts(simu_sce_ref1) <- log1p(counts(simu_sce_ref1))

saveRDS(simu_sce_ref1,"simu_sce_ref1.rds")

set.seed(123)
compare_figure_ref1 <- plot_reduceddim(ref_sce = example_sce1, 
                                  sce_list = list(simu_sce_ref1), 
                                  name_vec = c("Reference", "scDesign3"),
                                  assay_use = "logcounts", 
                                  if_plot = FALSE, 
                                  color_by = "celltype", 
                                  n_pc = 20)

compare_figure_ref1$ind <- rep(colData(example_sce1)[,"ind"],2)

saveRDS(compare_figure_ref1,"compare_figure_ref1.rds")

ggplot(compare_figure_ref1, aes_string(x = "UMAP1", y = "UMAP2",color = "celltype")) + 
  geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + scale_color_discrete(type = "viridis") +
  facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 

ggplot(compare_figure_ref1, aes_string(x = "UMAP1", y = "UMAP2",color = "ind")) + 
  geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + scale_color_discrete(type = "viridis") +
  facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom")

plot_ind(compare_figure_celltype_ind,example_sce1,simu_sce1,inds = c("1_1","2_2"))

plot_ind_celltype <- function(compare_figure_celltype_ind,example_sce1,simu_sce1,inds,celltypes,plot.type){
  compare_figure_celltype_ind <- compare_figure_celltype_ind[which(compare_figure_celltype_ind$ind%in%inds &
                                                                   compare_figure_celltype_ind$celltype%in%celltypes),]
  if(plot.type == "pca"){
    ggplot(compare_figure_celltype_ind, aes_string(x = "PC1", y = "PC2",color = "ind")) + 
      geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + 
      facet_wrap(ind~Method,nrow = length(inds)) + theme(aspect.ratio = 1, legend.position = "bottom") + 
      guides(color = guide_legend(override.aes = list(size = 2,alpha = 1)))
  }
  else if(plot.type == "umap"){
    ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "ind")) + 
      geom_point(alpha = 0.5, size = 1) +
      facet_wrap(ind~Method,nrow = length(inds)) + theme(aspect.ratio = 1, legend.position = "bottom") 
  }
}

plot_ind_celltype(compare_figure_celltype_ind = compare_figure_ref1,
             example_sce1 = example_sce1,
             simu_sce1 = simu_sce1,
             inds = unique(compare_figure_celltype_ind$ind)[1:5],
             celltypes = unique(compare_figure_celltype_ind$celltype)[1],
             plot = "umap")

plot_ind_celltype(compare_figure_celltype_ind = compare_figure_celltype_ind,
                  example_sce1 = example_sce1,
                  simu_sce1 = simu_sce1,
                  inds = unique(compare_figure_celltype_ind$ind)[1:5],
                  celltypes = unique(compare_figure_celltype_ind$celltype)[1],
                  plot = "umap")



