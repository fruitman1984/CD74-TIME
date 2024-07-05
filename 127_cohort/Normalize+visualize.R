library(tidyverse)
library(Seurat)
library(harmony)
library(ggraph)

future::plan("multisession",workers=1)
seurat.list <- SplitObject(Normal_harmony, split.by = "orig.ident")
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
#load('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74TIME/seuratlist_all.RData')
features <- SelectIntegrationFeatures(object.list = seurat.list,nfeatures = 3000)
Normal_harmony <- merge(seurat.list[[1]], y = seurat.list[2 : length(seurat.list)],merge.data = T)
VariableFeatures(Normal_harmony) <- features
options(future.globals.maxSize = 100000 * 1024^6)
future::plan("multisession",workers=8)
Normal_harmony <- ScaleData(Normal_harmony, features = features, verbose = FALSE,vars.to.regress = c("nCount_RNA","percent.mt"))
Normal_harmony <- RunPCA(object = Normal_harmony,  verbose = F)
Normal_harmony <- RunHarmony(Normal_harmony, group.by.vars=c("orig.ident"))
ElbowPlot(Normal_harmony,ndims = 50)
dim.use <- 1:50
Normal_harmony <- FindNeighbors(Normal_harmony, dims = dim.use, reduction = "harmony")
Normal_harmony <- FindClusters(Normal_harmony, resolution = 0.8)
Normal_harmony <- RunUMAP(Normal_harmony, dims = dim.use, reduction = "harmony")

DimPlot(Normal_harmony,label = T,pt.size = 1,raster = T)

Normal_harmony_sample <- subset(Normal_harmony,downsample=400)
future::plan('multisession',workers=8)
Idents(Normal_harmony_sample) <- 'seurat_clusters'
all.markers <- FindAllMarkers(Normal_harmony_sample, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file="Major_marker_clusterqc240325.xls",
            sep="\t",quote = F,row.names = F)

###### annotation ----
Major.idents <- c('Mono_Macro','CD8T','NK','GMP','GMP','GMP',
                  'CD4T','Eryth','EMP','CD4T','Eryth',
                  'B','EMP','Mono_Macro','Neutro','GMP',
                  'EMP','Eryth','Plasma','Eryth','Mono_Macro',
                  'EMP','Mono_Macro','DC','EMP','MK',
                  'T_prolif','Mono_Macro','Mono_Macro','MSC','DC',
                  'B','CD4T','B','GMP')
length(Major.idents)
names(Major.idents) <- levels(Normal_harmony)
Normal_harmony <- RenameIdents(Normal_harmony, Major.idents)
Normal_harmony@meta.data$'celltype'<- Idents(Normal_harmony)

Normal_harmony@meta.data$celltype <- factor(Normal_harmony@meta.data$celltype,levels = c('Mono_Macro','Neutro','DC','B','Plasma','CD4T','CD8T','NK','T_prolif','MK','GMP','EMP','MSC','Eryth'))
Idents(Normal_harmony) <- 'celltype'
cols_celltype <- c('Mono_Macro'='#ee827c','Neutro'='#38a1db','DC'='#ffd900','B'='#00a497','Plasma'='#98d98e','CD4T'='#83ccd2','CD8T'='#f08300','NK'='#cca6bf','T_prolif'='#a6a5c4','MK'='#bc763c','GMP'='#eaf4fc','EMP'='#d6e9ca','MSC'='#84a2d4','Eryth'='#ede4cd')
DimPlot(Normal_harmony,label = T,cols = cols_celltype,raster = T,pt.size = 1)
ggsave('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74TIME/Fig1/Umap_celltype.pdf',width = 5.2,height = 4)

