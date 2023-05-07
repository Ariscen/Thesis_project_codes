library(ggplot2)
library(SingleCellExperiment)
setwd("C:/Future_Dpan/sc-eQTL/data_backup/")

inds <- readRDS("inds.rds")
res <- readRDS("res.rds")
res_ref <- readRDS("res_ref.rds")
ind_res <- readRDS("ind_res.rds")
ind_res_ref <- readRDS("ind_res_ref.rds")

ind_res <- data.frame(do.call(rbind,ind_res))
ind_res$ind <- inds
ind_res$Method <- "scDesign-pop"

ind_res_ref <- data.frame(do.call(rbind,ind_res_ref))
ind_res_ref$ind <- inds
ind_res_ref$Method <- "scDesign3"

ind_res_all <- rbind(ind_res,ind_res_ref)
ind_res_all$Method <- factor(ind_res_all$Method,levels = c("scDesign3","scDesign-pop"))

# ind_res_all$ind <- as.character(ind_res_all$ind)
# ind_res_all <- rbind(ind_res_all,c(res,"All_cell","scDesign-pop"))
# ind_res_all <- rbind(ind_res_all,c(res_ref,"All_cell","scDesign3"))
# ind_res_all$PCA_lisi <- as.numeric(ind_res_all$PCA_lisi)
# ind_res_all$PCA_entropy <- as.numeric(ind_res_all$PCA_entropy)
# ind_res_all$PCA_cms <- as.numeric(ind_res_all$PCA_cms)
# ind_res_all$UMAP_lisi <- as.numeric(ind_res_all$UMAP_lisi)
# ind_res_all$UMAP_entropy <- as.numeric(ind_res_all$UMAP_entropy)
# ind_res_all$UMAP_cms <- as.numeric(ind_res_all$UMAP_cms)
# ind_res_all$Method <- factor(ind_res_all$Method,levels = c("scDesign3","scDesign-pop"))

ggplot(ind_res_all)+geom_boxplot(aes(Method,UMAP_lisi))+geom_jitter(aes(Method,UMAP_lisi))
ggplot(ind_res_all)+geom_boxplot(aes(Method,PCA_lisi))+geom_jitter(aes(Method,PCA_lisi))

shapiro.test(ind_res$PCA_lisi)
shapiro.test(ind_res_ref$PCA_lisi)
shapiro.test(ind_res$UMAP_lisi)
shapiro.test(ind_res_ref$UMAP_lisi)

t.test(ind_res$PCA_lisi,ind_res_ref$PCA_lisi)
wilcox.test(ind_res$PCA_lisi,ind_res_ref$PCA_lisi,paired = T)

t.test(ind_res$PCA_entropy,ind_res_ref$PCA_entropy)
wilcox.test(ind_res$PCA_entropy,ind_res_ref$PCA_entropy)

t.test(ind_res$PCA_cms,ind_res_ref$PCA_cms)
wilcox.test(ind_res$PCA_cms,ind_res_ref$PCA_cms)

t.test(ind_res$UMAP_lisi,ind_res_ref$UMAP_lisi)
wilcox.test(ind_res$UMAP_lisi,ind_res_ref$UMAP_lisi,paired = T)

t.test(ind_res$UMAP_entropy,ind_res_ref$UMAP_entropy)
wilcox.test(ind_res$UMAP_entropy,ind_res_ref$UMAP_entropy)

t.test(ind_res$UMAP_cms,ind_res_ref$UMAP_cms)
wilcox.test(ind_res$UMAP_cms,ind_res_ref$UMAP_cms)


compare_figure_combined <- readRDS("compare_figure_combined.rds")
example_sce1 <- readRDS("example_sce1.rds")
simu_sce_ref1 <- readRDS("simu_sce_ref1.rds")
simu_sce1 <- readRDS("simu_sce1.rds")

plot_all<-function(compare_figure,color_by){
  compare_figure <- compare_figure[which(compare_figure$Method!="scDesign3"),]
  ggplot(compare_figure, aes_string(x = "UMAP1", y = "UMAP2",color = color_by)) + 
    geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + scale_color_discrete(type = "viridis") +
    facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_all(compare_figure_combined,"celltype")

plot_all_ind<-function(compare_figure,color_by,inds){
  compare_figure <- compare_figure[which(compare_figure$Method!="scDesign3"),]
  compare_figure <- compare_figure[which(compare_figure$ind %in% inds),]
  ggplot(compare_figure, aes_string(x = "UMAP1", y = "UMAP2",color = color_by)) + 
    geom_point(alpha = 0.5, size = 1,aes_string(shape = NULL)) + scale_color_discrete(type = "viridis") +
    facet_wrap(Method~ind,ncol = length(inds)) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_all_ind(compare_figure_combined,"celltype",inds = inds[1:9])
plot_all_ind(compare_figure_combined,"celltype",inds = inds[10:17])

plot_feature <- function(compare_figure,example_sce1,simu_sce_ref1, simu_sce1,feature){
  compare_figure$feature <- c(logcounts(example_sce1)[feature,],
                              logcounts(simu_sce_ref1)[feature,],
                              logcounts(simu_sce1)[feature,])
  ggplot(compare_figure, aes_string(x = "UMAP1", y = "UMAP2",color = "feature")) + 
    geom_point(alpha = 0.5, size = 1) + scale_color_continuous(type = "viridis") +
    facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom") 
}

plot_feature(compare_figure_combined,example_sce1,simu_sce_ref1,simu_sce1,"LYZ")

snp_dataframe <- readRDS("snp_dataframe.rds")
eqtl_geno <- readRDS("eqtl_geno.rds")

plot_geno <- function(compare_figure,example_sce1,simu_sce_ref1, simu_sce1,feature,celltype,snp_dataframe){
  rsID <- as.character(snp_dataframe[which(snp_dataframe$gene_name==feature),"rsID"])
  #celltypes <- as.character(snp_dataframe[which(snp_dataframe$gene_name==feature),"cell_type"])
  tmp <- snp_dataframe[which(snp_dataframe$gene_name==feature),7:23]
  tmp1 <- tmp[compare_figure$ind]
  compare_figure$geno <- as.numeric(tmp1)
  
  # example_sce1 <- scater::logNormCounts(example_sce1)
  # simu_sce_ref1 <- scater::logNormCounts(simu_sce_ref1)
  # simu_sce1 <- scater::logNormCounts(simu_sce1)
  
  compare_figure$feature <- c(logcounts(example_sce1)[feature,],
                              logcounts(simu_sce_ref1)[feature,],
                              logcounts(simu_sce1)[feature,])
  compare_figure <- compare_figure[,c("Method","celltype","ind","geno","feature")]
  
  compare_figure_tmp <- dplyr::summarise(dplyr::group_by(compare_figure,Method,celltype,ind,geno),mean(feature))
  compare_figure_tmp <- compare_figure_tmp[which(compare_figure_tmp$celltype%in%celltype),]
  compare_figure_tmp <- data.frame(compare_figure_tmp)
  colnames(compare_figure_tmp)[5] <- "mean_expression"
  
  # model1 <- lm(mean_expression~geno,data = compare_figure_tmp[which(compare_figure_tmp$Method=="Reference"),])
  # print(summary(model1)$coefficient)
  # model2 <- lm(mean_expression~geno,data = compare_figure_tmp[which(compare_figure_tmp$Method=="scDesign3"),])
  # print(summary(model2)$coefficient)
  # model3 <- lm(mean_expression~geno,data = compare_figure_tmp[which(compare_figure_tmp$Method=="scDesign-pop"),])
  # print(summary(model3)$coefficient)
  
  ref <- compare_figure_tmp[which(compare_figure_tmp$Method=="Reference"),]
  print(cor.test(ref$mean_expression,ref$geno,method = "spearman"))
  scd3 <- compare_figure_tmp[which(compare_figure_tmp$Method=="scDesign3"),]
  print(cor.test(scd3$mean_expression,scd3$geno,method = "spearman"))
  scdp <- compare_figure_tmp[which(compare_figure_tmp$Method=="scDesign-pop"),]
  print(cor.test(scdp$mean_expression,scdp$geno,method = "spearman"))
  
  
  # ref <- as.numeric(compare_figure_tmp[which(compare_figure_tmp$Method=="Reference"),"mean_expression"])
  # scd3 <- as.numeric(compare_figure_tmp[which(compare_figure_tmp$Method=="scDesign3"),"mean_expression"])
  # scdp <- as.numeric(compare_figure_tmp[which(compare_figure_tmp$Method=="scDesign-pop"),"mean_expression"])
  # 
  # print(cor(ref,scd3,method = "pearson"))
  # print(cor(ref,scdp,method = "pearson"))
  
  ggplot(compare_figure_tmp)+geom_boxplot(aes(as.character(geno),mean_expression,group = geno))+
    ggtitle(label = rsID)+
    facet_wrap(celltype~Method,ncol = 3) + theme(aspect.ratio = 1, legend.position = "bottom",)+
    xlab(label = "Genotype")+ylab(label = "Expression level")
  
  
  # ggplot(compare_figure, aes_string(x = "UMAP1", y = "UMAP2",color = "geno")) + 
  #   geom_point(alpha = 0.5, size = 1) + scale_color_continuous(type = "viridis") +
  #   facet_wrap(~Method,nrow = 1) + theme(aspect.ratio = 1, legend.position = "bottom",title = rsID) 
}

unique(compare_figure_combined$celltype)

plot_geno(compare_figure_combined,
          example_sce1,
          simu_sce_ref1,
          simu_sce1,
          "BLK",
          "CD4 Naive/Central memory T cell",
          snp_dataframe)

plot_geno(compare_figure_combined,
          example_sce1,
          simu_sce_ref1,
          simu_sce1,
          "BLK",
          "Memory B Cell",
          snp_dataframe)

plot_geno(compare_figure_combined,
          example_sce1,
          simu_sce_ref1,
          simu_sce1,
          "LYZ",
          "Classic Monocyte",
          snp_dataframe)

