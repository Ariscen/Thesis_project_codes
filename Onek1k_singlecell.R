library(data.table)
library(R.utils)
setwd("~/projects/sceQTLsim/data/")

library(Seurat)
library(SeuratObject)
sc_data <- readRDS("~/projects/sceQTLsim/data/onek1k_analyzed/local.rds")
unique(sc_data@meta.data$cell_type_ontology_term_id)
unique(sc_data@meta.data$cell_type)
sort(unique(sc_data@meta.data$donor_id))

inds <- c("1_1","2_2","3_3","4_4","6_6","7_7","8_8","9_9","10_10","11_11","12_12","13_13",
          "15_15","16_16","17_17","18_18","19_19")
eqtl_geno <- eqtl_geno[,-which(!colnames(eqtl_geno)%in%c(inds,colnames(eqtl_geno)[1:6]))]

sc_data_17 <- subset(sc_data,donor_id%in%inds)
rm(sc_data)
table(sc_data_17@meta.data$cell_type)
DimPlot(sc_data_17,reduction="pca",group.by = "cell_type")

levels(sc_data_17@meta.data$cell_type)

celltypes <- c("natural killer cell","plasmacytoid dendritic cell","memory B cell","naive B cell",
               "gamma-delta T cell","regulatory T cell","transitional stage B cell",
               "naive thymus-derived CD4-positive, alpha-beta T cell","naive thymus-derived CD8-positive, alpha-beta T cell",
               "central memory CD4-positive, alpha-beta T cell","effector memory CD4-positive, alpha-beta T cell",
               "central memory CD8-positive, alpha-beta T cell","effector memory CD8-positive, alpha-beta T cell",
               "CD4-positive, alpha-beta cytotoxic T cell","CD16-negative, CD56-bright natural killer cell, human",
               "mucosal invariant T cell","plasmablast","conventional dendritic cell","CD14-positive monocyte",
               "CD14-low, CD16-positive monocyte")

sc_data_17_1 <- subset(sc_data_17,cell_type%in%celltypes)
unique(sc_data_17_1@meta.data$cell_type)

counts<-sc_data_17_1@assays$RNA@data

rownames(counts)<-gsub("_","-",as.character(sc_data_17_1@assays$RNA@meta.features$feature_name))
# sc_data_17_1@assays$RNA@meta.features$feature_id<-as.factor(rownames(sc_data_17_1))
# sc_data_17_1@assays$RNA@counts@Dimnames[[1]]<-gsub("_","-",as.character(sc_data_17_1@assays$RNA@meta.features$feature_name))
# sc_data_17_1@assays$RNA@data@Dimnames[[1]]<-gsub("_","-",as.character(sc_data_17_1@assays$RNA@meta.features$feature_name))

sc_data_17_2 <- CreateSeuratObject(counts)

sc_data_17_2 <- NormalizeData(sc_data_17_2, normalization.method = "LogNormalize", scale.factor = 10000)

sc_data_17_2 <- FindVariableFeatures(sc_data_17_2, selection.method = "vst", nfeatures = 2000)
sc_data_17_2 <- ScaleData(sc_data_17_2, features = rownames(sc_data_17_2))
sc_data_17_2 <- RunPCA(sc_data_17_2, features = VariableFeatures(object = sc_data_17_2))
DimPlot(sc_data_17_2, reduction = "pca")

ElbowPlot(sc_data_17_2,ndims = 30)
sc_data_17_2 <- FindNeighbors(sc_data_17_2, dims = 1:25)
sc_data_17_2 <- FindClusters(sc_data_17_2, resolution = 1)

sc_data_17_2 <- RunUMAP(sc_data_17_2, dims = 1:25)
DimPlot(sc_data_17_2, reduction = "umap",label = T)
FeaturePlot(sc_data_17_2,c("CD4","CD8A"),reduction = "umap")

sc_data_17_2@meta.data$original_celltype <- as.factor(as.character(sc_data_17_1@meta.data$cell_type))
sc_data_17_2@meta.data$donor_id <- as.factor(as.character(sc_data_17_1@meta.data$donor_id))
sc_data_17_2@meta.data$pool_number <- as.factor(as.character(sc_data_17_1@meta.data$pool_number))
sc_data_17_2@meta.data$age <- as.factor(as.character(sc_data_17_1@meta.data$age))

table(sc_data_17_2@meta.data$original_celltype[which(sc_data_17_2@active.ident=="9")])

VlnPlot(sc_data_17_2,c("GZMA","GZMB","XCL1","XCL2"))
VlnPlot(sc_data_17_2,c("CD14","LYZ","FCGR3A","FCGR3B"))
VlnPlot(sc_data_17_2,c("TCL1A","TNFRSF13B","CD27","TNFRSF17"))
VlnPlot(sc_data_17_2,c("CST3","FCER1A","SERPINF1"))
VlnPlot(sc_data_17_2,c("S100B","IL7R","GZMK"))
VlnPlot(sc_data_17_2,c("GNLY","NKG7"))
VlnPlot(sc_data_17_2,c("LTB","CCR7","PASK"))
VlnPlot(sc_data_17_2,c("SOX4","ID2","SELL","CCR7"))
FeaturePlot(sc_data_17_2,c("SOX4","ID2","SELL"))
FeaturePlot(sc_data_17_2,c("KLRB1","CD8A","CD4"))
VlnPlot(sc_data_17_2,c("GZMK","TNFSF13B","KLRB1"))
VlnPlot(sc_data_17_2,c("CCR7","SELL","LRRN3"))
DimPlot(sc_data_17_2, reduction = "umap",label = T)
# [1] "Memory B Cell"                   "CD4 Naive/Central memory T cell" "CD8 Effector memory"            
# [4] "CD8 Naive/Central memory T cell" "Natural Killer Cell"             "Non-classic Monocyte"           
# [7] "CD4 Effector memory/TEMRA"       "Classic Monocyte"                "CD8 S100B T cell"               
# [10] "Naïve/Immature B Cell"           "Dendritic Cell"                  "Natural Killer Recruiting Cell" 
# [13] "CD4 SOX4 T cell"                 "Plasma Cell"

new.cluster.ids <- c(
  "CD4 Naive/Central memory T cell",
  "CD4 Effector memory/TEMRA",
  "CD4 Naive/Central memory T cell",
  "CD8 Effector memory",
  "Natural Killer Cell",
  "CD8 Effector memory",
  "Naïve/Immature B Cell",
  "Memory B Cell",
  "CD4 SOX4 T cell" ,
  "CD8 Naive/Central memory T cell",
  "Classic Monocyte",
  "Natural Killer Recruiting Cell",
  "Non-classic Monocyte",
  "Plasma Cell",
  "Dendritic Cell",
  "CD8 S100B T cell",
  "Dendritic Cell"
)
names(new.cluster.ids) <- levels(sc_data_17_2)
sc_data_17_2 <- RenameIdents(sc_data_17_2, new.cluster.ids)
sc_data_17_2@active.ident<-factor(sc_data_17_2@active.ident,
                                  levels = sort(unique(eqtl_geno$cell_type)))
DimPlot(sc_data_17_2, reduction = "umap",label = T)

saveRDS(sc_data_17_2,"sc_data_17_2.rds")

unique(sc_data_17_2@active.ident)%in%unique(eqtl_geno$cell_type)
unique(eqtl_geno$cell_type)
# sc_data@assays$RNA@counts<-sc_data@assays$RNA@data
# 
# sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
# sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
# 
# Inf%in%sc_data@assays$RNA@data@x
# 
# top10 <- head(VariableFeatures(sc_data), 10)
# sc_data@assays$RNA@meta.features$feature_name[which(rownames(sc_data)%in%top10)]
# plot1 <- VariableFeaturePlot(sc_data)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

# sc_data <- ScaleData(sc_data, features = rownames(sc_data))
# sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))
# sc_data <- GetAssayData(sc_data,slot = "data")
# sc_data <- CreateSeuratObject(counts = sc_data, project = "onek1k", 
#                            min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
# Visualize QC metrics as a violin plot
# VlnPlot(sc_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# try predicted.celltype.l1 to start 
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# # remotes::install_github("mojaveazure/seurat-disk")
# library(SeuratDisk)
# library(ggplot2)
# library(patchwork)
# 
# reference <- LoadH5Seurat("~/projects/sceQTLsim/data/pbmc_multimodal.h5seurat")
# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", 
#         label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# 
# #devtools::install_github('satijalab/seurat-data')
# library(SeuratData)
# #InstallData('pbmc3k')
# sc_data@assays$RNA@counts <- sc_data@assays$RNA@data
# 
# sc_data <- SCTransform(sc_data, verbose = TRUE)
# 
# anchors <- FindTransferAnchors(
#   reference = reference,
#   query = sc_data,
#   normalization.method = "SCT",
#   reference.reduction = "spca",
#   dims = 1:50
# )
# 
# sc_data <- MapQuery(
#   anchorset = anchors,
#   query = sc_data,
#   reference = reference,
#   refdata = list(
#     celltype.l1 = "celltype.l1",
#     celltype.l2 = "celltype.l2",
#     predicted_ADT = "ADT"
#   ),
#   reference.reduction = "spca", 
#   reduction.model = "wnn.umap"
# )
# 
# p1 = DimPlot(sc_data, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# p2 = DimPlot(sc_data, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
# p1 + p2

# try reclustering from the start
# as.character(unique(sc_data@meta.data$predicted.celltype.l2))
# type_selected<-setdiff(as.character(unique(sc_data@meta.data$predicted.celltype.l2)),c("Doublet","Eryth","Platelet"))
# sc_data<-subset(sc_data,predicted.celltype.l2 %in% type_selected)
# as.character(unique(sc_data@meta.data$predicted.celltype.l2))
# sc_data@assays$RNA@counts<-sc_data@assays$RNA@data
# 
# sc_data<-NormalizeData(sc_data,)

sc_data_ind1<-subset(sc_data,donor_id == "1_1")
#sc_data_ind2<-subset(sc_data,donor_id == "2_2")
sc_data_ind2<-subset(sc_data,donor_id == "691_692")

sc_data@meta.data$donor_id[1]
table(sc_data@meta.data$cell_type)
table(sc_data@meta.data$predicted.celltype.l2)

table(sc_data_ind2@meta.data$predicted.celltype.l2)
table(sc_data_ind2@meta.data$cell_type)

sc_data@assays$RNA@data@Dimnames[[1]]<-as.character(sc_data@assays$RNA@meta.features$feature_name)

FeaturePlot(sc_data,"XCL1",reduction = "umap",raster = FALSE)
DimPlot(sc_data,reduction = "umap",group.by = "predicted.celltype.l2",label = T)
DimPlot(sc_data,reduction = "umap",group.by = "cell_type",raster = FALSE,label = T)

cell_type_number<-read.csv("~/projects/sceQTLsim/cell_type_numbers.csv",header = T)
sum(cell_type_number[which(cell_type_number$Cell.Type=="CD4 NC"),"Number.of.Cells"])

tmps<-c()
for( i in 1:length(unique(cell_type_number$Individual))){
  tmp <- unique(cell_type_number[which(cell_type_number$Individual==i),"Total.Number.of.Cells.per.Person"])
  tmps <- c(tmps,tmp)
}






