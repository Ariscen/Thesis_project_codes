library(data.table)
library(R.utils)
library(Seurat)
library(pbapply)
setwd("~/projects/sceQTLsim/data/")
sc_data_17_2 <- readRDS("sc_data_17_2.rds")
DimPlot(sc_data_17_2,label = T)
eqtl_geno <- readRDS("eqtl_geno.rds")
gene_overlap <- unique(eqtl_geno$gene_name)

# gene
sc <- GetAssayData(sc_data_17_2,slot = "counts")[gene_overlap,]
dim(sc)

# celltype
celltype <- factor(as.character(sc_data_17_2@active.ident),
                   levels = sort(unique(as.character(sc_data_17_2@active.ident))))
levels(celltype)

# individual id
ind <- factor(sc_data_17_2@meta.data$donor_id,
              levels = c("1_1","2_2","3_3","4_4","6_6","7_7","8_8","9_9","10_10","11_11","12_12","13_13",
                         "15_15","16_16","17_17","18_18","19_19"))
ind

# merge into dataframe
data <- data.frame(ind=ind, celltype=celltype)
rownames(data) <- colnames(sc)
data$gene <- t(as.matrix(sc))

saveRDS(data,"data.rds")

# snps
gene_snps <- pblapply(gene_overlap,function(g){
  eqtl_tmp <- eqtl_geno[which(eqtl_geno$gene_name==g),]
  snps_tmp <- unique(eqtl_tmp$rsID)
  snp_geno <- lapply(snps_tmp,function(snp){
    geno_tmp1 <- eqtl_tmp[which(eqtl_tmp$rsID==snp)[1],levels(ind)]
    names(geno_tmp1) <- levels(ind)
    geno_tmp2 <- geno_tmp1[ind]
    return(geno_tmp2)
  })
  names(snp_geno) <- snps_tmp
  return(snp_geno)
})
names(gene_snps) <- gene_overlap

saveRDS(gene_snps,"gene_snps.rds")




