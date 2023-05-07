library(data.table)
library(R.utils)
setwd("~/projects/sceQTLsim/data/")

library(vcfR)
d_chr1<-read.vcfR("~/projects/sceQTLsim/data/onek1k/vcfs/filter_vcf_r08/chr1.dose.filtered.R2_0.8.vcf.gz")

library(magrittr)
dirpath <- "/home/chrisycd/proj/scEQTLsim"
inputdir <- sprintf("%s/data/OneK1K/analysis", dirpath)

# specify col type for ID else read in error
cols_spec <- readr::cols(ID = "c")

# read in sample ids
eqtl_geno <- sprintf("%s/all_eqtl_geno.tsv", inputdir) %>% 
  readr::read_delim(., delim = "\t", col_types = cols_spec, 
                    col_names = TRUE)
inds <- c("1_1","2_2","3_3","4_4","6_6","7_7","8_8","9_9","10_10","11_11","12_12","13_13",
          "15_15","16_16","17_17","18_18","19_19")
eqtl_geno <- eqtl_geno[,-which(!colnames(eqtl_geno)%in%c(inds,colnames(eqtl_geno)[1:6]))]

# get overlap with 2000 HVG
library(Seurat)
sc_data_17_2 <- readRDS("sc_data_17_2.rds")
gene_overlap <- intersect(unique(eqtl_geno$gene_name),VariableFeatures(sc_data_17_2))

eqtl_geno <- eqtl_geno[which(eqtl_geno$gene_name%in%gene_overlap),] 

saveRDS(eqtl_geno,"eqtl_geno.rds")
