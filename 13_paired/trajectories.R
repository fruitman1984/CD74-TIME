######### slingshot --------
setwd('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid/slingshot')
# BiocManager::install("slingshot")
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(ggrastr)
library(tradeSeq)
library(RColorBrewer)

Myeloid_Mono <- subset(Myeloid,idents = c('cDC_CD1C','pDC_JCHAIN'),invert=T)
Myeloid_Mono@meta.data$celltype2 <- factor(Myeloid_Mono@meta.data$celltype2,levels = c('Mono_MPO','Mono_STAB1','Mono_PRKCD','Mono_MS4A4E',
                                                                                       'Mono_GRN','Mono_NOTCH2','Mono_MKI67','Macro_CLEC6A','Macro_LY86','Macro_LGALS2','Macro_LY6E','Macro_CD74'))

col_Myeloid <- c('Mono_MPO'='#cee4ae','Mono_STAB1'='#84b9cb','Mono_PRKCD'='#5b7e91','Mono_MS4A4E'='#b68d4c',
                 'Mono_GRN'='#5383c3','Mono_NOTCH2'='#69821b','Mono_MKI67'='#c7b370',
                 'Macro_CLEC6A'='#65318e','Macro_LY86'='#bb5535','Macro_LGALS2'='#efab93',
                 'Macro_LY6E'='#74325c','Macro_CD74'='#d7003a')

col_celltype <- col_Myeloid[Myeloid_Mono$celltype2]
names(col_celltype) <- Myeloid_Mono$celltype2
sce <- as.SingleCellExperiment(Myeloid_Mono)
Myeloid_Mono_slingshot <- slingshot(sce, clusterLabels = "celltype2", reducedDim = "SCVI_UMAP",start.clus = c('Mono_MPO'),stretch = 0)
summary(Myeloid_Mono_slingshot$slingPseudotime_1)

Myeloid_Mono@tools[['slingshot']] = SlingshotDataSet(Myeloid_Mono_slingshot)
pseudotime = slingPseudotime(Myeloid_Mono_slingshot)
nc <- 3
pt <- slingPseudotime(Myeloid_Mono_slingshot)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)

for (i in nms) {
  pseudotime_sub <- pseudotime[colnames(Myeloid_Mono),i]
  Myeloid_Mono <- AddMetaData(object = Myeloid_Mono,
                              metadata = pseudotime_sub,
                              col.name = i
  )
}

umap <- reducedDims(Myeloid_Mono_slingshot)$SCVI_UMAP %>% as.data.frame() %>% cbind(celltype = Myeloid_Mono_slingshot$celltype2) %>% cbind(pt)
curves <- slingCurves(Myeloid_Mono_slingshot,as.df = TRUE)
ggplot(umap,aes(x= scviumap_1 , y = scviumap_2)) +
  labs(title="Slingshot")+
  geom_point_rast(aes(color = celltype), size = 0.3)+
  scale_color_manual(values = col_Myeloid)+
  geom_path(data = curves,aes(group = Lineage), linewidth = 0.5)+
  theme( plot.title = element_text(hjust = 0.5),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         panel.background = element_rect(fill = 'white'),
         plot.background=element_rect(fill="white"),
         axis.line.x = element_line(color="black", linewidth = 0.5),
         axis.text.x = element_text(colour = 'black',size = 10),
         axis.text.y = element_text(colour = 'black',size = 10),
         axis.line.y = element_line(color="black", linewidth = 0.5),
         legend.title = element_blank(),
         legend.key.size=unit(0.5,'cm'),
         legend.key=element_rect(fill='white')) +
  guides(color = guide_legend(override.aes = list(size=2))) 
ggsave('slingshot_umap_paired_mono.pdf',width = 4.5,height = 4)

average <- apply(pt, 1, function(x) mean(x, na.rm = TRUE))
umap$peseudotime <- average

ggplot(data = umap,aes(x= scviumap_1 , y = scviumap_2)) +
  geom_point_rast(aes(color = peseudotime), size = 0.3)+
  scale_colour_viridis_c(na.value="#D3D3D3", option = "C")+
  #scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name ="Spectral")))(100))+
  geom_path(data = curves,aes(group = Lineage), linewidth = 0.5)+
  theme( plot.title = element_text(hjust = 0.5),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         panel.background = element_rect(fill = 'white'),
         plot.background=element_rect(fill="white"),
         axis.line.x = element_line(color="black", linewidth = 0.5),
         axis.text.x = element_text(colour = 'black',size = 10),
         axis.text.y = element_text(colour = 'black',size = 10),
         axis.line.y = element_line(color="black", linewidth = 0.5),
         legend.title = element_blank(),
         legend.key.size=unit(0.5,'cm'),
         legend.key=element_rect(fill='white'))
ggsave('slingshot_pseudotime_umap.pdf',width = 4,height = 3.5)

####### monocle3 --------
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
setwd('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid/monocle3')
data <- GetAssayData(Myeloid_Mono, assay = 'RNA', slot = 'counts')
cell_metadata <- Myeloid_Mono@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Myeloid_Mono, reduction = "scvi_umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

col_Myeloid <- c('Mono_MPO'='#cee4ae','Mono_STAB1'='#84b9cb','Mono_PRKCD'='#5b7e91','Mono_MS4A4E'='#b68d4c',
                 'Mono_GRN'='#5383c3','Mono_NOTCH2'='#69821b','Mono_MKI67'='#c7b370',
                 'Macro_CLEC6A'='#65318e','Macro_LY86'='#bb5535','Macro_LGALS2'='#efab93',
                 'Macro_LY6E'='#74325c','Macro_CD74'='#d7003a')

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "celltype2",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,cell_size=0.8)+ scale_color_manual(values=col_Myeloid)

cds <- order_cells(cds)
plot_ras <- plot_cells(cds,
                       color_cells_by = "pseudotime",
                       label_cell_groups=FALSE,
                       label_leaves=F,
                       label_branch_points=F,
                       graph_label_size=1.5,show_trajectory_graph = F,
                       group_label_size=4,cell_size=0.5,rasterize = T)+scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
plot_ras <- ggrastr::rasterise(plot_ras,layer='Point',dpi =300)
plot_ras

plot_ras <- plot_cells(cds,
                       color_cells_by = "pseudotime",
                       label_cell_groups=FALSE,
                       label_leaves=T,
                       label_branch_points=T,
                       graph_label_size=1.5,show_trajectory_graph = T,
                       group_label_size=4,cell_size=0.5,rasterize = T)+scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
plot_ras <- ggrastr::rasterise(plot_ras,layer='Point',dpi =300)
plot_ras


###### cytotrace ------
setwd('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid/cytotrace')
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
library(CytoTRACE2)

Myeloid_Mono_exp <- GetAssayData(Myeloid_Mono,assay = 'RNA',slot = 'data')
cytotrace2_result <- cytotrace2(Myeloid_Mono_exp,species = 'human')

Myeloid_Mono <- AddMetaData(Myeloid_Mono,metadata = cytotrace2_result)
save(cytotrace2_result,file = 'cytotrace2_result.RData')

FeaturePlot(Myeloid_Mono,features = 'CytoTRACE2_Score',raster = T,pt.size = 3.5)+scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name ="Spectral")))(100))
ggsave('CytoTRACE2_Score_Myeloid.pdf',width = 4.6,height = 4)

FeaturePlot(Myeloid_Mono,features = 'CytoTRACE2_Relative',raster = T,pt.size = 3.5)+scale_colour_gradientn(colours = colorRampPalette(brewer.pal(n = 10, name ="Spectral"))(100))
ggsave('CytoTRACE2_Score_Myeloid_rev.pdf',width = 4.6,height = 4)