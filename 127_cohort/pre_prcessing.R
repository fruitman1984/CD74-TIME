library(Seurat)
library(tidyverse)
library(ggraph)
library(harmony)
setwd('/GPUFS/sysu_zhaom_1/QG/300/85')

############# Batch1 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/old/Batch1")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/old/Batch1/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch1 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch1@meta.data$Batch <- 'Batch1'

save(seurat_Batch1,file = 'seurat_Batch1.RData')
rm(list=ls())
gc()

############# Batch2 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/old/Batch2")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/old/Batch2/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch2 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch2@meta.data$Batch <- 'Batch2'

save(seurat_Batch2,file = 'seurat_Batch2.RData')
rm(list=ls())
gc()

############# Batch3 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/old/Batch3")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/old/Batch3/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch3 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch3@meta.data$Batch <- 'Batch3'

save(seurat_Batch3,file = 'seurat_Batch3.RData')
rm(list=ls())
gc()

############# Batch4 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix_64252")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix_64252/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch4 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch4@meta.data$Batch <- 'Batch4'

save(seurat_Batch4,file = 'seurat_Batch4.RData')
rm(list=ls())
gc()

############# Batch5 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix2")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix2/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch5 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch5@meta.data$Batch <- 'Batch5'

save(seurat_Batch5,file = 'seurat_Batch5.RData')
rm(list=ls())
gc()


################
load('/GPUFS/sysu_zhaom_1/QG/300/85/seurat_Batch1.RData')
load('/GPUFS/sysu_zhaom_1/QG/300/85/seurat_Batch2.RData')
load('/GPUFS/sysu_zhaom_1/QG/300/85/seurat_Batch3.RData')
load('/GPUFS/sysu_zhaom_1/QG/300/85/seurat_Batch4.RData')
load('/GPUFS/sysu_zhaom_1/QG/300/85/seurat_Batch5.RData')

seurat_data <- merge(seurat_Batch5,list(seurat_Batch4,seurat_Batch3,seurat_Batch2,seurat_Batch1),merge.data = T)

seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
VlnPlot(seurat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0,group.by = 'orig.ident')
gene.freq <- do.call("cbind", tapply(seurat_data@meta.data$nFeature_RNA,seurat_data@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(seurat_data@meta.data$nCount_RNA,seurat_data@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(seurat_data@meta.data$percent.mt,seurat_data@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))

write.csv(freq.combine,file='freq_combine_qc.csv')

######## scDblFinder
load('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_data_merge.RData')
library(scDblFinder)
library(BiocParallel) 
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seurat_data_merge) 
sce <- scDblFinder(sce, samples="orig.ident") 

table(sce$scDblFinder.class)
scDbl_data <- sce@colData[,c("scDblFinder.class","scDblFinder.score","scDblFinder.weighted","scDblFinder.cxds_score")] %>% as.data.frame()
seurat_data_merge <- AddMetaData(seurat_data_merge,scDbl_data)


######## decontX 
library(decontX)
counts <- seurat_data_merge@assays$RNA@counts
decontX_results <- decontX(counts) 
save(decontX_results,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/Batch_decontX_results.RData')
seurat_data_merge$Contamination <- decontX_results$contamination

seurat_data_merge <- subset(seurat_data_merge,nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 1000 & percent.mt < 10 & Contamination < 0.2)

seurat_data_merge <- subset(seurat_data_merge,nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 20)
save(seurat_data_merge,file = '/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_data_merge.RData')

############# Batch6 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix240301")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix240301/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch6 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch6@meta.data$Batch <- 'Batch6'

save(seurat_Batch6,file = 'seurat_Batch6.RData')
rm(list=ls())
gc()

load('/GPUFS/sysu_zhaom_1/QG/300/127/seurat_Batch6.RData')
library(scDblFinder)
library(BiocParallel) 
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seurat_Batch6) 
sce <- scDblFinder(sce, samples="orig.ident") 
table(sce$scDblFinder.class)
scDbl_data <- sce@colData[,c("scDblFinder.class","scDblFinder.score","scDblFinder.weighted","scDblFinder.cxds_score")] %>% as.data.frame()
seurat_Batch6 <- AddMetaData(seurat_Batch6,scDbl_data)
save(seurat_Batch6,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch6.RData')


library(decontX)
counts <- seurat_Batch6@assays$RNA@counts
decontX_results <- decontX(counts) 
save(decontX_results,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/Batch6_decontX_results.RData')
seurat_Batch6$Contamination <- decontX_results$contamination
save(seurat_Batch6,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch6.RData')

############# Batch7 --------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix240701")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix240701/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch7 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch7@meta.data$Batch <- 'Batch7'

save(seurat_Batch7,file = '/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch7.RData')
rm(list=ls())
gc()

load('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch7.RData')
library(scDblFinder)
library(BiocParallel) 
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seurat_Batch7) 
sce <- scDblFinder(sce, samples="orig.ident") 

table(sce$scDblFinder.class)
scDbl_data <- sce@colData[,c("scDblFinder.class","scDblFinder.score","scDblFinder.weighted","scDblFinder.cxds_score")] %>% as.data.frame()
seurat_Batch7 <- AddMetaData(seurat_Batch7,scDbl_data)
save(seurat_Batch7,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch7.RData')

pdf()
FeaturePlot(seurat_Batch7,features = "scDblFinder.score")
dev.off()

seurat_Batch7 <- AddMetaData(seurat_Batch7,scDbl_data)

library(decontX)
counts <- seurat_Batch7@assays$RNA@counts
decontX_results <- decontX(counts) 
save(decontX_results,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/Batch7_decontX_results.RData')
seurat_Batch7$Contamination <- decontX_results$contamination
save(seurat_Batch7,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch7.RData')

############ Batch 8 ---------
samples <- dir(path="/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix240628")
seurat_list <- list()
for (i in samples){
  seurat_list[[i]] <- assign(i, Read10X(data.dir = paste0("/GPUFS/sysu_zhaom_1/QG/300/matrix/matrix240628/", i)))
  colnames(seurat_list[[i]]) <- paste0(i,'-',colnames(seurat_list[[i]]))
}

for (i in seq_along(seurat_list)){
  seurat_list[[i]] <- CreateSeuratObject(
    seurat_list[[i]],
    project = "SeuratProject", 
    min.cells = 10,
    min.features = 200,
    names.field = 1,
    names.delim = "-")
}

seurat_Batch8 <- merge(seurat_list[[1]], y = seurat_list[2 : length(seurat_list)],merge.data = T)

seurat_Batch8@meta.data$Batch <- 'Batch8'

save(seurat_Batch8,file = '/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch8.RData')
rm(list=ls())
gc()

load('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch8.RData')
library(scDblFinder)
library(BiocParallel) 
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seurat_Batch8) 
sce <- scDblFinder(sce, samples="orig.ident") 

table(sce$scDblFinder.class)
scDbl_data <- sce@colData[,c("scDblFinder.class","scDblFinder.score","scDblFinder.weighted","scDblFinder.cxds_score")] %>% as.data.frame()
seurat_Batch8 <- AddMetaData(seurat_Batch8,scDbl_data)

pdf()
FeaturePlot(seurat_Batch8,features = "scDblFinder.score")
dev.off()

seurat_Batch8 <- AddMetaData(seurat_Batch8,scDbl_data)

library(decontX)
counts <- seurat_Batch8@assays$RNA@counts
decontX_results <- decontX(counts)
save(decontX_results,file='/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/Batch8_decontX_results.RData')
seurat_Batch8$Contamination <- decontX_results$contamination
save(seurat_Batch8,file = '/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_Batch8.RData')

######## QC -----
setwd('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData')
seurat_data_merge2 <- merge(seurat_Batch6,y=list(seurat_Batch7,seurat_Batch8),merge.data = T)
seurat_data_merge2[["percent.mt"]] <- PercentageFeatureSet(seurat_data_merge2, pattern = "^MT-")
save(seurat_data_merge2,file = '/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74_240701/RData/seurat_data_merge2.RData')


VlnPlot(seurat_data_merge2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0,group.by = 'orig.ident')
gene.freq <- do.call("cbind", tapply(seurat_data_merge2@meta.data$nFeature_RNA,seurat_data_merge2@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(seurat_data_merge2@meta.data$nCount_RNA,seurat_data_merge2@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(seurat_data_merge2@meta.data$percent.mt,seurat_data_merge2@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
write.csv(freq.combine,file='freq_combine_qc2.csv')

VlnPlot(seurat_data_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0,group.by = 'orig.ident')
gene.freq <- do.call("cbind", tapply(seurat_data_merge@meta.data$nFeature_RNA,seurat_data_merge@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(seurat_data_merge@meta.data$nCount_RNA,seurat_data_merge@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(seurat_data_merge@meta.data$percent.mt,seurat_data_merge@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
write.csv(freq.combine,file='freq_combine_qc.csv')

seurat_data_merge2 <- subset(seurat_data_merge2,nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 1000 & percent.mt < 10 & scDblFinder.class == 'singlet' & Contamination < 20)
seurat_data_merge <- subset(seurat_data_merge,nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 1000 & percent.mt < 10 & scDblFinder.class == 'singlet' & Contamination < 20)
