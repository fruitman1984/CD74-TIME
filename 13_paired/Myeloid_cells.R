library(Seurat)
library(reticulate)
library(tidyverse)

setwd('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid')
load('/media/user/sdg/home/wqh/data_copy/IME/1/Myeliod.RData')

adata = sc$read_h5ad('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myeloid_h.h5ad')
scvi_embed <- adata$obsm['X_scVI']
scvi_embed <- as.matrix(scvi_embed)
rownames(scvi_embed) <- py_to_r(adata$obs_names$to_frame())[[1]]
Myeloid[["scvi"]] <- CreateDimReducObject(embeddings = scvi_embed, key = "scvi_", assay = DefaultAssay(Myeloid))
Myeloid <- FindNeighbors(Myeloid, dims = 1:30, reduction = "scvi")
Myeloid <- FindClusters(Myeloid, resolution =1)
Myeloid <- RunUMAP(Myeloid, dims = 1:30, reduction = "scvi", n.components = 2,reduction.name = 'scvi_umap')
DimPlot(Myeloid, reduction = "scvi_umap", pt.size =2,raster = T,label = T)

Myeloid_sample <- subset(Myeloid,downsample=500)
future::plan('multisession',workers=16)
options(future.globals.maxSize = 100000 * 1024^6)
DefaultAssay(Myeloid_sample) <- 'RNA'
Idents(Myeloid_sample) <- 'seurat_clusters'
all.markers <- FindAllMarkers(Myeloid_sample, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file="Myeloid_scvi_RNA.xls",
            sep="\t",quote = F,row.names = F)

Major.idents <- c('Mono_STAB1','Mono_PRKCD','Macro_CLEC6A','cDC_CD1C','Mono_MKI67','Mono_MS4A4E',
                  'Macro_LY86','Mono_GRN','Macro_CD74','Mono_NOTCH2','Macro_LGALS2',
                  'Macro_LY6E','Mono_MPO','pDC_JCHAIN','Mono_PRKCD')
length(Major.idents)
names(Major.idents) <- levels(Myeloid)
Myeloid <- RenameIdents(Myeloid, Major.idents)
Myeloid@meta.data$'celltype2'<- Idents(Myeloid)

col_Myeloid <- c('Mono_MPO'='#cee4ae','Mono_STAB1'='#84b9cb','Mono_PRKCD'='#5b7e91','Mono_MS4A4E'='#b68d4c',
                 'Mono_GRN'='#5383c3','Mono_NOTCH2'='#69821b','Mono_MKI67'='#c7b370',
                 'Macro_CLEC6A'='#65318e','Macro_LY86'='#bb5535','Macro_LGALS2'='#efab93',
                 'Macro_LY6E'='#74325c','Macro_CD74'='#d7003a','cDC_CD1C'='#316745','pDC_JCHAIN'='#9e8b8e')

Myeloid@meta.data$celltype2 <- factor(Myeloid@meta.data$celltype2,levels = c('Mono_MPO','Mono_STAB1','Mono_PRKCD','Mono_MS4A4E',
                                                                             'Mono_GRN','Mono_NOTCH2','Mono_MKI67','Macro_CLEC6A','Macro_LY86','Macro_LGALS2','Macro_LY6E','Macro_CD74',
                                                                             'cDC_CD1C','pDC_JCHAIN'))
Idents(Myeloid) <- Myeloid@meta.data$celltype2
DimPlot(Myeloid,label = T,raster = T,pt.size = 2.5,cols = col_Myeloid,label.box = T,label.size = 3)
ggsave('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid/Umap_paired_Myeloid_celltype.pdf',width = 6,height = 4.5)

cell.prop<-as.data.frame(prop.table(table(Myeloid$celltype2, Myeloid$status),margin = 2))
colnames(cell.prop)<-c("celltype","group","proportion")
cell.prop$proportion <- cell.prop$proportion*100
ggplot(cell.prop,aes(group,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab('Proportion (%)')+
  theme_bw()+
  theme(axis.ticks.length=unit(0.25,'cm'),axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(colour = 'black'))+
  guides(fill=guide_legend(title=NULL))+ scale_fill_manual(values=col_Myeloid)
ggsave("Proportion_myeloid_celltype_status.pdf",width = 3,height = 4)

patient1 <- subset(Myeloid@meta.data,status == 'RE')
cell.prop<-table(patient1$celltype2)
cell.prop1 <- cell.prop %>% as.data.frame()
colnames(cell.prop1) <- c('celltype','RE')
patient1 <- subset(Myeloid@meta.data,status == 'DX')
cell.prop<-table(patient1$celltype2)
cell.prop2 <- cell.prop %>% as.data.frame()
cell.prop1$DX <- cell.prop2$Freq
cell.prop1.1 <- cell.prop1[1:4,]
colSums(cell.prop1.1[2:3])
cell.prop1.2 <- cell.prop1[5,]
cell.prop1.3 <- data.frame(celltype='Others',RE=4097,DX=5865)
cell.prop1.x <- rbind(cell.prop1.2,cell.prop1.3)
result <- fisher.test(cell.prop1.x[, c("DX", "RE")])
result
write.csv(cell.prop1.x,file = 'Myeloid_chisquared.csv')

######## scProgram ------
library(scProgram)
source('/media/user/sdg/home/qiuguo/10X/TME/scProgram.R')

FeatureMatrix = GetFeatures(obj = Myeloid, group.by = "celltype2", genenumber = 50, pct_exp = 0.1, mode = "fast")
FeatureMatrix <- FeatureMatrix[,levels(Myeloid$celltype2)]
HeatFeatures(obj = Myeloid, features = FeatureMatrix, group.by = "celltype2", 
             show_rownames = F, show_colnames = T, cols = c("white","white", "white", "#52A85F"))
# 4*5
unique(Idents(Myeloid))
Myeloid_Mono <- subset(Myeloid,idents = c('cDC_CD1C','pDC_JCHAIN'),invert=T)
Myeloid_Mono@meta.data$celltype2 <- factor(Myeloid_Mono@meta.data$celltype2,levels = c('Mono_MPO','Mono_STAB1','Mono_PRKCD','Mono_MS4A4E',
                                                                                       'Mono_GRN','Mono_NOTCH2','Mono_MKI67','Macro_CLEC6A','Macro_LY86','Macro_LGALS2','Macro_LY6E','Macro_CD74'))
FeatureMatrix = GetFeatures(obj = Myeloid_Mono, group.by = "celltype2", genenumber = 50, pct_exp = 0.1, mode = "fast")
FeatureMatrix <- FeatureMatrix[,levels(Myeloid_Mono$celltype2)]
HeatFeatures(obj = Myeloid_Mono, features = FeatureMatrix, group.by = "celltype2", 
             show_rownames = F, show_colnames = T, cols = c("white","white", "white", "#52A85F"))

GetProgram2 <- function(features = FeatureMatrix,
                        genesets = c("KEGG", "GO"), # 可以是任意组合的基因集
                        pvalue_cutoff = 0.05,
                        cols = c("#F47E5D", "#CA3D74", "#7F2880", "#463873"),
                        plot_term_number = 3) {
  
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
  gmtfiles <- list(
    KEGG = system.file("data", "h.all.v7.2.symbols.c2.cp.kegg.v7.2.symbols.gmt", package = "scProgram"),
    GO = '/media/user/sdg/home/qiuguo/datasets/GSEA/c5.go.bp.v2023.1.Hs.symbols.gmt',
    REACTOME = '/media/user/sdg/home/qiuguo/datasets/GSEA/c2.cp.reactome.v2023.1.Hs.symbols.gmt'
  )
  
  c5 <- do.call(rbind, lapply(genesets, function(gs) {
    read.gmt(gmtfiles[[gs]])
  }))
  
  ego_comb <- data.frame()
  cluster.names <- colnames(FeatureMatrix)
  
  for (cluster.name in cluster.names) {
    gene_list <- FeatureMatrix[, cluster.name]
    ego <- enricher(gene_list, TERM2GENE = c5)
    ego_result <- ego@result
    ego_result$GeneRatio.num <- ego_result$Count / length(gene_list)
    ego_result <- ego_result[order(ego_result$GeneRatio.num, decreasing = TRUE), ]
    ego_result <- subset(ego_result, ego_result$pvalue < pvalue_cutoff)
    ego_result$group <- cluster.name
    ego_comb <- rbind(ego_comb, ego_result)
  }
  
  ego_comb_sub <- ego_comb %>%
    group_by(group) %>%
    slice_max(order_by = -log10(pvalue), n = plot_term_number) #
  pal <- colorRampPalette(cols)(100)
  ego_comb_sub$logp <- -log10(ego_comb_sub$pvalue)
  ego_comb_sub$Description <- factor(ego_comb_sub$Description, levels = rev(unique(ego_comb_sub$Description)))
  ego_comb_sub$group <- factor(ego_comb_sub$group, levels = rev(unique(ego_comb_sub$group)))
  
  return(ego_comb_sub)
}

FeatureMatrix = GetFeatures(obj = Myeloid_Mono, group.by = "celltype2", genenumber = 100, pct_exp = 0.1, mode = "fast")
FeatureMatrix <- FeatureMatrix[,levels(Myeloid_Mono$celltype2)]

ego_comb_sub <- GetProgram2(features = FeatureMatrix, pvalue_cutoff = 0.05,
                            cols = c("#F47E5D", "#CA3D74", "#7F2880", "#463873"), plot_term_number =50)
write.csv(ego_comb_sub,file = '/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid/ego_com_sub.csv',row.names = F)
ego_comb_sub <- read.csv('/media/user/sdg/home/qiuguo/10X/TME/NEW/Myleoid/Myeloid_pathway_selected.csv',header = T)

ego_comb_sub$Description <- factor(ego_comb_sub$Description, levels = rev(unique(ego_comb_sub$Description)))
ego_comb_sub$group <- factor(ego_comb_sub$group, levels = levels(Myeloid_Mono$celltype2))
ggplot(ego_comb_sub, aes(x = group, y = Description)) +
  geom_point(aes(size = Count, color = logp)) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  scale_color_gradientn(colours = pal <- colorRampPalette(rev(c("#F47E5D", "#CA3D74", "#463873")))(100)) + #"#7F2880",
  ylab(NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1))+ #aspect.ratio=1,
  ggtitle("")


######## Ro/e -------
library(tidyverse)
library(readr)
library(circlize)
library(ComplexHeatmap)
ROIE <- function(crosstab){
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

crosstab_Myeloid <- table(Myeloid@meta.data[['celltype2']], Myeloid@meta.data[['status']])
roie_matrix <- ROIE(crosstab_Myeloid)

breakpoints <- c(min(roie_matrix), 2, max(roie_matrix))
colors <- c("white", "#E67E22", "#E67E22")
my_col <- colorRamp2(breakpoints, colors)

Heatmap(roie_matrix, name = "Ro/e", col = my_col,cluster_columns = F,cluster_rows = F, border = NA, cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", roie_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
})

####### Macrophage -----
unique(Idents(Myeloid))
Myeloid_Macro <- subset(Myeloid,idents = c('Macro_CLEC6A','Macro_LY86','Macro_LGALS2','Macro_LY6E','Macro_CD74'))

Genes <- c('ITGAM','CD68','CD14','FCGR1A','CD80','CD86','IL23A','CCR7','KYNU','IDO1','CD40','IRF1','IRF5','IL1A','IL1B','IL6','CCL5','CXCL9','CXCL10','CXCL11','TNF',
           'CD163','MRC1','FCER2','FCGR2A','IRF4','FN1','MSR1','CD276','TNFSF8','TNFSF12','WNT7B','CLEC7A','MMP9','MMP14','MMP19','TGFB1','TGFB2','TGFB3','CTSA','CTSB','CTSD','EGF','VEGFA','VEGFB',
           'LYVE1','CCL4','CCL18','CCL20','CCL22','IL1RN','IL1R2','IL4R','IL10','CSF1R','CD274','PDCD1LG2','ARG1','ARG2','CD200R1')
Myeloid_Macro_ave <- AverageExpression(Myeloid_Macro,features = Genes,group.by = 'celltype2',return.seurat = F,assays = 'RNA',slot = 'data')[[1]] %>% as.data.frame()

anno_row <- data.frame(anno=c(rep('M0',3),rep('M1',18),rep('M2',39)),row.names = Genes)
pheatmap(Myeloid_Macro_ave,gaps_row = c(3,21),cluster_rows = F,cluster_cols = F,scale = 'row',annotation_row = anno_row,border_color = NA)

library(reticulate)
library(sceasy)
use_condaenv('scvi-env')
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
adata <- convertFormat(Myeloid_Macro, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)
adata$write_h5ad(filename = 'Myeloid_Macro_rna.h5ad')

adata = sc$read_h5ad('Myeloid_Macro_2.h5ad')
scvi_embed <- adata$obsm['X_scVI']
scvi_embed <- as.matrix(scvi_embed)
rownames(scvi_embed) <- py_to_r(adata$obs_names$to_frame())[[1]]
Myeloid_Macro[["scvi"]] <- CreateDimReducObject(embeddings = scvi_embed, key = "scvi_", assay = DefaultAssay(Myeloid_Macro))

Myeloid_Macro <- FindNeighbors(Myeloid_Macro, dims = 1:30, reduction = "scvi")
Myeloid_Macro <- FindClusters(Myeloid_Macro, resolution =1.5)
Myeloid_Macro <- RunUMAP(Myeloid_Macro, dims = 1:30, reduction = "scvi", n.components = 2,reduction.name = 'scvi_umap')
DimPlot(Myeloid_Macro, reduction = "scvi_umap", pt.size =2,raster = T,label = T,split.by = 'status')
DotPlot(Myeloid_Macro,features = c('HAVCR2','CXCL13','PDCD1','LAYN','TOX','IFNG','MIR155HG','TNFRSF9','ITGAE','CD8A','CD4'))

Myeloid_Macro@meta.data$celltype2 <- factor(Myeloid_Macro@meta.data$celltype2,levels = c('Macro_CLEC6A','Macro_LY86','Macro_LGALS2','Macro_LY6E','Macro_CD74'))
col_Macro <- c('Macro_CLEC6A'='#65318e','Macro_LY86'='#bb5535','Macro_LGALS2'='#efab93','Macro_LY6E'='#74325c','Macro_CD74'='#d7003a')
DimPlot(Myeloid_Macro, reduction = "scvi_umap", pt.size =8,raster = T,label = F,raster.dpi = c(1500,1500),split.by = 'status',cols = col_Macro,group.by = 'celltype2')
ggsave('Umap_Macro_paired_split.pdf',width = 7.5,height = 4)

cell.prop<-as.data.frame(prop.table(table(Myeloid_Macro$celltype2, Myeloid_Macro$patient),margin = 2))
colnames(cell.prop)<-c("celltype","group","proportion")
ggplot(cell.prop,aes(group,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.25,'cm'),axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(colour = 'black'))+
  guides(fill=guide_legend(title=NULL))+ scale_fill_manual(values=col_Macro)
ggsave("Proportion_macro_celltype_patient.pdf",width = 8,height = 4)

cell.prop<-as.data.frame(prop.table(table(Myeloid_Macro$celltype2, Myeloid_Macro$status),margin = 2))
colnames(cell.prop)<-c("celltype","group","proportion")
cell.prop$proportion <- cell.prop$proportion*100
ggplot(cell.prop,aes(group,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ylab('Proportion (%)')+
  theme_bw()+
  theme(axis.ticks.length=unit(0.25,'cm'),axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(colour = 'black'))+
  guides(fill=guide_legend(title=NULL))+ scale_fill_manual(values=col_Macro)
ggsave("Proportion_macro_celltype_status.pdf",width = 3,height = 4)

patient1 <- subset(Myeloid_Macro@meta.data,status == 'RE')
cell.prop<-table(patient1$celltype2)
cell.prop1 <- cell.prop %>% as.data.frame()
colnames(cell.prop1) <- c('celltype','RE')
patient1 <- subset(Myeloid_Macro@meta.data,status == 'DX')
cell.prop<-table(patient1$celltype2)
cell.prop2 <- cell.prop %>% as.data.frame()
cell.prop1$DX <- cell.prop2$Freq
cell.prop1.1 <- cell.prop1[1:4,]
colSums(cell.prop1.1[2:3])
cell.prop1.2 <- cell.prop1[5,]
cell.prop1.3 <- data.frame(celltype='Others',RE=3333,DX=2560)
cell.prop1.x <- rbind(cell.prop1.2,cell.prop1.3)
result <- fisher.test(cell.prop1.x[, c("DX", "RE")])
result
write.csv(cell.prop1.x,file = 'Macro_chisquared.csv')