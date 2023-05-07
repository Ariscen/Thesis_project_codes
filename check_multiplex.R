library(stringr)
setwd("/home/ycen/projects/sceQTLsim/data/onek1k/GenotypeSamples/")
ls <- list.files(pattern = "*.gz")

# pre-check
tmps <- c()
for (l in ls){
  tmp <- read.table(gzfile(l))
  print(as.vector(as.matrix(tmp)))
  tmps <- c(tmps,as.vector(as.matrix(tmp)))
}

table(tmps)[which(table(tmps) > 1)]
length(unique(tmps))
tmps_uni <- unique(tmps)

# collect sample id
samp <- c()
for (i in 1:length(ls)){
  samp <- c(samp,str_split(ls,"_")[[i]][4])
}
samp

# create a matrix to record individuals in each sample
ind_samp <- matrix(data = 0, nrow = length(tmps_uni), ncol = length(samp))
rownames(ind_samp) <- tmps_uni
colnames(ind_samp) <- samp
for (i in 1:length(ls)){
  vec <- as.vector(as.matrix(read.table(gzfile(ls[i]))))
  for (v in vec){
    ind_samp[v,i] = ind_samp[v,i] + 1 
  }
}
summary(ind_samp)
saveRDS(ind_samp,"/home/ycen/projects/sceQTLsim/onek1k_individuals_samples_mat.rds")

# check multiplex
rowSums(ind_samp)[rowSums(ind_samp)>1]
ind_samp <- readRDS("/home/ycen/projects/sceQTLsim/onek1k_individuals_samples_mat.rds")
which(ind_samp==2)
ind_samp[which(rownames(ind_samp)=="643_644"),]
# sample19：869_870 and 870_871; sample38: 870_871; sample76: 564_656
s19 <- read.table(gzfile("/home/ycen/projects/sceQTLsim/data/onek1k/GenotypeSamples/GSM5899891_OneK1K_scRNA_Sample19_GenotypeSamples.txt.gz"))
