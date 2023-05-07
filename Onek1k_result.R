library(scDesign3)
library(ggplot2)
library(SingleCellExperiment)
theme_set(theme_bw())

setwd("C:/Future_Dpan/sc-eQTL/data_backup/")

data <- readRDS("data.rds")
snp_matrix <- readRDS("snp_matrix.rds")
marginal_list <- readRDS("marginal_list.rds")

dim(counts(example_sce))
example_sce1 <-SingleCellExperiment(assays = list(counts = t(data$gene)),colData = data[,1:2])
example_sce1

example_sce1 <- readRDS("example_sce1.rds")

dat <- data.frame(celltype = data$celltype, ind = data$ind, corr_group = as.numeric(data$ind))
Onek1k_data <- list(count_mat = data$gene, dat = dat, geno = snp_matrix, NewCovariate = NULL)

Onek1k_data <- readRDS("Onek1k_data.rds")

Onek1k_copula <- fit_copula(
  sce = example_sce1,
  assay_use = "counts",
  marginal_list = marginal_list,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 2,
  new_covariate = NULL,
  input_data = Onek1k_data$dat
)

Onek1k_copula <- readRDS("Onek1k_copula.rds")

Onek1k_para <- extract_para(
  sce = example_sce1,
  marginal_list = marginal_list,
  n_cores = 1,
  family_use = "nb",
  new_covariate = NULL
)

Onek1k_newcount <- simu_new(
  sce = example_sce1,
  mean_mat = Onek1k_para$mean_mat,
  sigma_mat = Onek1k_para$sigma_mat,
  zero_mat = Onek1k_para$zero_mat,
  quantile_mat = NULL,
  copula_list = Onek1k_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = Onek1k_data$dat,
  new_covariate = Onek1k_data$new_covariate,
  important_feature = Onek1k_copula$important_feature
)

logcounts(example_sce1) <- log1p(counts(example_sce1))
simu_sce1 <- example_sce1
counts(simu_sce1) <- Onek1k_newcount
logcounts(simu_sce1) <- log1p(counts(simu_sce1))
saveRDS(simu_sce1,"simu_sce1.rds")

set.seed(123)
compare_figure_celltype <- plot_reduceddim(ref_sce = example_sce1, 
                                           sce_list = list(simu_sce1), 
                                           name_vec = c("Reference", "scDesign3-pop"),
                                           assay_use = "logcounts", 
                                           if_plot = FALSE, 
                                           color_by = "celltype", 
                                           n_pc = 20)
mean(compare_figure_celltype$celltype[1:20146] == Onek1k_data$dat$celltype)


compare_figure_ind <- plot_reduceddim(ref_sce = example_sce1, 
                                      sce_list = list(simu_sce1), 
                                      name_vec = c("Reference", "scDesign3-pop"),
                                      assay_use = "logcounts", 
                                      if_plot = FALSE, 
                                      color_by = "ind", 
                                      n_pc = 20)
mean(compare_figure_ind$ind[1:20146] == Onek1k_data$dat$ind)


compare_figure_celltype_ind <- compare_figure_celltype
compare_figure_celltype_ind$ind <- compare_figure_ind$ind
compare_figure_celltype_ind <- compare_figure_celltype_ind[,c("Method","celltype","ind",
                                                              colnames(compare_figure_celltype[3:24]))]

saveRDS(compare_figure_celltype_ind, "compare_figure_celltype_ind.rds")
compare_figure_celltype_ind <- readRDS("compare_figure_celltype_ind.rds")

ggplot(compare_figure_celltype_ind, aes_string(x = "PC1", y = "PC2",color = "celltype")) + 
  geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + 
  facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size = 2,alpha = 1)))

ggplot(compare_figure_celltype_ind, aes_string(x = "PC1", y = "PC2",color = "ind")) + 
  geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + 
  facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size = 2,alpha = 1)))

ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "celltype")) + 
  geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + 
  facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 

ggplot(compare_figure_celltype_ind, aes_string(x = "UMAP1", y = "UMAP2",color = "ind")) + 
  geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + 
  facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 
