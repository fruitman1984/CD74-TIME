library(Seurat)
library(reticulate)
library(tidyverse)
library(harmony)

load('seurat_paired_raw.RData')

seurat_paired[["percent.mt"]] <- PercentageFeatureSet(seurat_paired, pattern = "^MT-")
VlnPlot(seurat_paired, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0,group.by = 'patient')
gene.freq <- do.call("cbind", tapply(seurat_paired@meta.data$nFeature_RNA,seurat_paired@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(seurat_paired@meta.data$nCount_RNA,seurat_paired@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(seurat_paired@meta.data$percent.mt,seurat_paired@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))

write.csv(freq.combine,file='freq_combine_qc.csv')

seurat_paired <- subset(seurat_paired,nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA > 1000 & nCount_RNA < 35000 & percent.mt < 20)
save(seurat_paired,file = 'seurat_paired.RData')

library(harmony)
future::plan("multisession",workers=1)
seurat.list <- SplitObject(seurat_paired, split.by = "orig.ident")
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seurat.list,nfeatures = 3000)
seurat_paired <- merge(seurat.list[[1]], y = seurat.list[2 : length(seurat.list)],merge.data = T)
VariableFeatures(seurat_paired) <- features
options(future.globals.maxSize = 100000 * 1024^6)
future::plan("multisession",workers=3)
seurat_paired <- ScaleData(seurat_paired, features = features, verbose = FALSE,vars.to.regress = c("nCount_RNA","percent.mt"))
seurat_paired <- RunPCA(object = seurat_paired,  verbose = F)
seurat_paired <- RunHarmony(seurat_paired, group.by.vars=c("orig.ident"))
ElbowPlot(seurat_paired,ndims = 50)
dim.use <- 1:30
seurat_paired <- FindNeighbors(seurat_paired, dims = dim.use, reduction = "harmony")
seurat_paired <- FindClusters(seurat_paired, resolution = 0.6)
seurat_paired <- RunUMAP(seurat_paired, dims = dim.use, reduction = "harmony")

DimPlot(seurat_paired,raster=T,pt.size=2)

##### predict-----
setwd('/media/user/sdg/home/qiuguo/10X/TME/NEW/predicted_based/predicted')

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(randomForest)


## load classifier
load("/media/user/sdg/home/qiuguo/10X/TME/NEW/predicted_based/predicted/RF1.RData", verbose = TRUE)   # as generated using the RFtrain script
load("/media/user/sdg/home/qiuguo/10X/TME/NEW/predicted_based/predicted/RF2.RData", verbose = TRUE)

# X <- log(E+1)
X <- as.matrix(GetAssayData(seurat_paired,assay = 'RNA',slot = 'data'))
## predict, within seconds
# RF1

subX <- X[intersect(rownames(X),rownames(rf.all.inner$importance)),]
unmatched_genes <- rownames(rf.all.inner$importance)[!(rownames(rf.all.inner$importance) %in% rownames(X))]
unmatched_df <- data.frame(matrix(0, nrow = length(unmatched_genes), ncol = ncol(subX)),row.names = unmatched_genes)
colnames(unmatched_df) <- colnames(subX)
subX <- rbind(subX,unmatched_df)

X.predict.all <- predict(rf.all.inner, newdata = t(subX), type="prob")

X.predict.all.max <- apply(X.predict.all, 1, max)
X.predict.all.class <- colnames(X.predict.all)[apply(X.predict.all, 1, which.max)]
table(X.predict.all.class)

# RF2

subX <- X[intersect(rownames(X),rownames(rf.tumor.inner$importance)),]
unmatched_genes <- rownames(rf.tumor.inner$importance)[!(rownames(rf.tumor.inner$importance) %in% rownames(X))]
unmatched_df <- data.frame(matrix(0, nrow = length(unmatched_genes), ncol = ncol(subX)),row.names = unmatched_genes)
colnames(unmatched_df) <- colnames(subX)
subX <- rbind(subX,unmatched_df)

X.predict.tumor <- predict(rf.tumor.inner, newdata = t(subX), type="prob")

X.predict.tumor.max <- apply(X.predict.tumor, 1, max)
X.predict.tumor.class <- colnames(X.predict.tumor)[apply(X.predict.tumor, 1, which.max)]
table(X.predict.tumor.class)

## write output
d.out <- data.frame(RF1.class=X.predict.all.class, RF1.score=X.predict.all.max, RF1=X.predict.all, 
                    RF2.class=X.predict.tumor.class, RF2.score=X.predict.tumor.max, RF2=X.predict.tumor)

write.table(d.out, file="predict.txt", quote = FALSE, sep = "\t")

save(list=c("X.predict.all", "X.predict.tumor"), file="predict.RData")

##### predict outs 
d.out <- read.table('/media/user/sdg/home/qiuguo/10X/TME/NEW/predicted_based/predicted/predict.txt',header = T,row.names = 1,sep='\t')
seurat_paired <- AddMetaData(seurat_paired,metadata = d.out)
seurat_paired@meta.data$CellType_RF <- ifelse(grepl("^Tumor\\.", seurat_paired@meta.data$RF2.class), "AML", "Normal")
DimPlot(seurat_paired, group.by = 'CellType_RF',cols = c('firebrick2','skyblue2'),raster = T,raster.dpi = c(1500,1500),pt.size = 3)
ggsave('Umap_paired_harmony_CellType_RF.pdf',width = 5,height = 4)

seurat_paired@meta.data$CellType1 <- ifelse(seurat_paired@meta.data$CellType_RF == 'AML','AML',seurat_paired@meta.data$RF2.class)
DimPlot(seurat_paired,label = T,raster = T,pt.size = 1,group.by = 'CellType1')
DimPlot(seurat_paired,label = T,group.by = 'CellType1',cols = col_celltype,raster = T,pt.size = 6,raster.dpi = c(2000,2000))
ggsave('Umap_harmony_CellType1.pdf',width = 5,height = 4)

DimPlot(seurat_paired,label = T,raster = T,pt.size = 1,group.by = 'seurat_clusters')
ggsave('Umap_harmony_clusters.pdf',width = 5.5,height = 4)

Idents(seurat_paired) <- 'CellType1'

######### TIME ---------
TIME_paired <- subset(seurat_paired,idents= 'AML',invert=T)

TIME_paired <- FindVariableFeatures(TIME_paired)
options(future.globals.maxSize = 100000 * 1024^6)
future::plan("multisession",workers=8)
TIME_paired <- ScaleData(TIME_paired, verbose = FALSE,vars.to.regress = c("nCount_RNA","percent.mt"))
TIME_paired <- RunPCA(object = TIME_paired,  verbose = F)
TIME_paired <- RunHarmony(TIME_paired, group.by.vars=c("orig.ident"))
ElbowPlot(TIME_paired,ndims = 50)
dim.use <- 1:30
TIME_paired <- FindNeighbors(TIME_paired, dims = dim.use, reduction = "harmony")
TIME_paired <- FindClusters(TIME_paired, resolution = c(0.8))
TIME_paired <- RunUMAP(TIME_paired, dims = dim.use, reduction = "harmony")
save(TIME_paired,file = 'TIME_paired_new.RData')
DimPlot(TIME_paired,label = T,pt.size = 2,raster = T,split.by = 'status')

TIME_paired_sample <- subset(TIME_paired,downsample=500)
future::plan('multisession',workers=16)
options(future.globals.maxSize = 100000 * 1024^6)
DefaultAssay(TIME_paired_sample) <- 'RNA'
Idents(TIME_paired_sample) <- 'seurat_clusters'
all.markers <- FindAllMarkers(TIME_paired_sample, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file="TIME_paired_RNA.xls",
            sep="\t",quote = F,row.names = F)