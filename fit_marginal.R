library(data.table)
library(R.utils)
library(Seurat)
library(pbapply)
library(pbmcapply)
library(scDesign3)
setwd("~/projects/sceQTLsim/data/")

data <- readRDS("data.rds")
gene_snps <- readRDS("gene_snps.rds")
eqtl_geno <- readRDS("eqtl_geno.rds")

data$ind
data$celltype
#gene <- "VIM"

genes <- colnames(data$gene)

snp_list <- pbmclapply(genes,function(gene){
  #print(gene)
  data1 <- data
  data1$gene <- data1$gene[,gene]
  if(length(gene_snps[[gene]])==1){
    gene_snps1 <- as.numeric(as.vector(gene_snps[[gene]][[1]]))
    return(gene_snps1)
  }else{
    gene_snps1 <- t(do.call(rbind,gene_snps[[gene]]))
    for(x in 1:length(gene_snps1[1,])){
      #print(length(unique(gene_snps1[,x])))
      if(length(unique(gene_snps1[,x]))==3){
        break
      }
    }
    gene_snps1 <- gene_snps1[,x]
    return(gene_snps1)
  }

},mc.cores = 2)


snp_matrix<-do.call(rbind,snp_list)
rownames(snp_matrix)<-colnames(data$gene)
snp_matrix <- t(snp_matrix)
saveRDS(snp_matrix,"snp_matrix.rds")


marginal_list <- pbmclapply(genes,function(gene){
  #print(gene)
  data1 <- data
  data1$gene <- data1$gene[,gene]
  if(length(gene_snps[[gene]])==1){
    gene_snps1 <- as.numeric(as.vector(gene_snps[[gene]][[1]]))
    data1 <- cbind(data1,gene_snps1)
    colnames(data1)[4] <- "genotype"
  }else{
    gene_snps1 <- t(do.call(rbind,gene_snps[[gene]]))
    for(x in 1:length(gene_snps1[1,])){
      #print(length(unique(gene_snps1[,x])))
      if(length(unique(gene_snps1[,x]))==3){
        break
      }
    }
    gene_snps1 <- gene_snps1[,x]
    data1 <- cbind(data1,gene_snps1)
    colnames(data1)[4] <- "genotype"
  }
  
  if(length(unique(data1$genotype))==3){
    mgcv.fit <- tryCatch({
      res <- mgcv::gam(formula = gene~celltype+
                         s(genotype, k = 3, bs = "cr")+
                         s(ind,bs="re"),
                       data = data1, family = "nb")
      res
    }, error = function(error) {
      message(paste0(gene, " gam fit fails!"))
      NULL
    }, silent = FALSE)
    
    return(mgcv.fit)
  }else{
    mgcv.fit <- tryCatch({
      res <- mgcv::gam(formula = gene~celltype+
                         #s(genotype, k = 2, bs = "cr")+
                         s(ind,bs="re"),
                       data = data1, family = "nb")
      res
    }, error = function(error) {
      message(paste0(gene, " gam fit fails!"))
      NULL
    }, silent = FALSE)
    
    return(mgcv.fit)
  }
},mc.cores = 2)

# summary(mgcv.fit)
# gam.vcomp(mgcv.fit)
# AIC(mgcv.fit)
# 
# par(mfrow = c(2,2))
# gam.check(mgcv.fit)
# plot(mgcv.fit)

saveRDS(marginal_list,"marginal_list.rds")


rm(list = ls())
