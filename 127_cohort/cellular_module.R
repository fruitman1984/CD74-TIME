setwd('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74TIME/module')

aaa <- prop.table(table(Normal_TIME$orig.ident,Normal_TIME$celltype2),margin = 1) %>% as.data.frame()
bbb <- reshape2::dcast(data = aaa,formula = Var1~Var2,fill = 'Freq')
bbb <- tibble::column_to_rownames(bbb,var='Var1')
bbb[,c(1:ncol(bbb))] <- as.numeric(unlist(bbb[,c(1:ncol(bbb))]))

df_list <- list()
data_list <- c(Normal_Myeloid,Normal_NK,Normal_Bcell,Normal_Tcell)
for (i in c(1:4)) {
  aaa <- prop.table(table(data_list[[i]]$orig.ident,data_list[[i]]$celltype2),margin = 1) %>% as.data.frame()
  bbb <- reshape2::dcast(data = aaa,formula = Var1~Var2,fill = 'Freq')
  bbb <- tibble::column_to_rownames(bbb,var='Var1')
  bbb[,c(1:ncol(bbb))] <- as.numeric(unlist(bbb[,c(1:ncol(bbb))]))
  df_list[[i]] <- bbb
}

for (i in c(1:4)) {
  df_list[[i]] <- df_list[[i]][unique(Normal_TIME$orig.ident),]
  df_list[[i]][is.na(df_list[[i]])] <- 0
}

combined_df <- bind_cols(df_list)
rownames(combined_df) <- unique(Normal_TIME$orig.ident)


combined_df_clinical <- combined_df
survival_metadata <- read.csv('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74TIME/module/class3/survival_metadata.csv',header = T,row.names = 1)
survival_metadata <- survival_metadata[rownames(combined_df_clinical),]
for (i in colnames(combined_df_clinical)) {
  combined_df_clinical[,paste0(i)] <- ifelse(combined_df_clinical[,i] > median(combined_df_clinical[,i]),1,0)
}
survival_metadata <- cbind(combined_df_clinical,survival_metadata)
survival_metadata$RFS <- ifelse(survival_metadata$event == 'Survival' & survival_metadata$REL == '0', 0,1)
survival_metadata$tRFS <- survival_metadata$tREL

M <- cor(combined_df,method = 'pearson')
pdf(paste0('Corrplot_eachtype_','pearson','_','ward.D2','.pdf'),width = 12,height = 12)
cor_dat <- corrplot(M, method = 'color', hclust.method='ward.D2',addrect = 7,order = 'hclust',col =colorRampPalette(colors = c('blue',"skyblue","white","firebrick3",'firebrick3'))(100),tl.cex=0.8,pch.col = 'black',tl.col = 'black')
dev.off()
pdf(paste0('Forrest_corrplot_','pearson','_','ward.D2','.pdf'),width = 8,height = 13)
forrest_plot(survival_metadata,ncol = 71,event = 'RFS',time = 'tRFS',order = colnames(cor_dat$corr))
dev.off()

forrest_order <- c('Plasma_IgM_IgA_IgG','CD8T_MKI67','CD8T_TYMS','Mac_CD74_MT1X','CD4T_CD52','CD8T_HAVCR2',
                   'cDC_CD1C','Plasma_IgA_IgG_SDC1','Mono_MPO','Mac_IL1B','Plasma_IgA_IgG_CD27','Mono_IFI30','Plasma_IgA_IgG_MKI67','Mac_TLR4',
                   'Plasma_IgA_IgG_IGLC1','B_IgM_IgD_FOXP1','NK_CREM','CD4T_SF1','Mono_CCNL1','NK_IKZF2','NK_NKTR','CD8T_FGFBP2',
                   'Mono_PLBD1', 'NK_CD56','Mono_VCAN','NK_CCL3','Mono_JUN','NK_IFNG','CD8T_CD226', 'B_CD24_IER2','B_IgM_IgD_CD69','CD8T_CD160','CD8T_CXCR4',
                   'B_CD34','NK_CD16','B_IgM_IgD_FCRL1','Plasma_IgM_IgA_CD27','B_IgM_IgD_CD24_CD1c','NK_PDCD4','Mac_LY6E','NK_HLADRA','CD8T_CCL5','pDC_LILRA4','pDC_JCHAIN', 'CD4T_ANXA1','NK_IFI6','CD8T_IFIT3', 
                   'CD4T_LTB','gdT_TRGC1_TRDC','Mac_CD74_CD49d','B_CD24_TXNIP','MAIT_SLC4A10','B_CD11c_HLADR','NK_GZMH', 'CD4T_FOXP3','Mono_S100A10','CD8T_LYAR','B_CD11b_ANXA2','CD8T_CD16',
                   'Mac_CD163','CD8T_LEF1','B_IgM_IgD_TCL1A','CD4T_MAL','NK_PRF1','Mono_CYP11B1','Neu_CD16b','NK_ADGRG1','Neu_CD66b','NK_ITGB1','NK_MT1X','Mono_APOC1')

pdf(paste0('Forrest2_corrplot_','pearson','_','ward.D2','.pdf'),width = 8,height = 13)
forrest_plot(survival_metadata,ncol = 71,event = 'RFS',time = 'tRFS',order = forrest_order)
dev.off()

Normal_TIME$module <- ifelse(Normal_TIME$celltype2 %in% c('CD4T_CD52','Mac_CD74_MT1X','Plasma_IgM_IgA_IgG','CD8T_HAVCR2','CD8T_MKI67','CD8T_TYMS'),'CM1',
                             ifelse(Normal_TIME$celltype2 %in% c('Mac_TLR4','Mac_IL1B','cDC_CD1C','Plasma_IgA_IgG_SDC1','Plasma_IgA_IgG_CD27',
                                                                 'Plasma_IgA_IgG_MKI67','Mono_IFI30','Mono_MPO'),'CM2',
                                    ifelse(Normal_TIME$celltype2 %in% c('NK_CREM','B_IgM_IgD_FOXP1','Mono_CCNL1','NK_NKTR','CD4T_SF1','CD8T_FGFBP2','CD8T_CD226', 'NK_IKZF2','Plasma_IgA_IgG_IGLC1'),'CM3',
                                           ifelse(Normal_TIME$celltype2 %in% c('CD8T_CXCR4','CD8T_CD160','NK_CD56','B_CD24_IER2','NK_CCL3','NK_IFNG', 'B_IgM_IgD_CD69','Mono_JUN','Mono_VCAN','Mono_PLBD1'),'CM4',
                                                  ifelse(Normal_TIME$celltype2 %in% c('Mac_LY6E', 'B_IgM_IgD_FCRL1','NK_IFI6','CD8T_IFIT3','B_CD34','pDC_LILRA4',
                                                                                      'pDC_JCHAIN','B_IgM_IgD_CD24_CD1c', 'Plasma_IgM_IgA_CD27','NK_CD16','NK_HLADRA','CD4T_ANXA1','NK_PDCD4','CD8T_CCL5'),'CM5',
                                                         ifelse(Normal_TIME$celltype2 %in% c('CD4T_FOXP3','B_CD11c_HLADR','CD4T_LTB', 'NK_GZMH','gdT_TRGC1_TRDC','CD8T_CD16','Mono_S100A10','Mac_CD74_CD49d','MAIT_SLC4A10','B_CD24_TXNIP','CD8T_LYAR','B_CD11b_ANXA2'),'CM6',
                                                                ifelse(Normal_TIME$celltype2 %in% c('NK_MT1X','Mono_APOC1','Mac_CD163','Neu_CD66b','Mono_CYP11B1','Neu_CD16b',
                                                                                                    'CD4T_MAL','CD8T_LEF1','NK_PRF1','B_IgM_IgD_TCL1A','NK_ITGB1','NK_ADGRG1'),'CM7','Others')))))))


aaa <- prop.table(table(Normal_TIME$orig.ident,Normal_TIME$module),margin = 1) %>% as.data.frame()
bbb <- reshape2::dcast(data = aaa,formula = Var1~Var2,fill = 'Freq')
bbb <- tibble::column_to_rownames(bbb,var='Var1')
bbb[,c(1:ncol(bbb))] <- as.numeric(unlist(bbb[,c(1:ncol(bbb))]))
survival_metadata <- read.csv('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74TIME/module/survival_metadata.csv',header = T,row.names = 1)

survival_metadata <- survival_metadata[rownames(bbb),]
# for (i in colnames(combined_df_clinical)) {
#   combined_df_clinical[,paste0(i)] <- ifelse(combined_df_clinical[,i] > median(combined_df_clinical[,i]),1,0)
# }
survival_metadata <- cbind(bbb,survival_metadata)
survival_metadata$RFS <- ifelse(survival_metadata$event == 'Survival' & survival_metadata$REL == '0', 0,1)
survival_metadata$tRFS <- survival_metadata$tREL

survival_metadata$CM1group <- ifelse(survival_metadata$CM1 > median(survival_metadata$CM1),1,0)
survival_metadata$CM2group <- ifelse(survival_metadata$CM2 > median(survival_metadata$CM2),1,0)
survival_metadata$CM3group <- ifelse(survival_metadata$CM3 > median(survival_metadata$CM3),1,0)
survival_metadata$CM4group <- ifelse(survival_metadata$CM4 > median(survival_metadata$CM4),1,0)
survival_metadata$CM5group <- ifelse(survival_metadata$CM5 > median(survival_metadata$CM5),1,0)
survival_metadata$CM6group <- ifelse(survival_metadata$CM6 > median(survival_metadata$CM6),1,0)
survival_metadata$CM7group <- ifelse(survival_metadata$CM7 > median(survival_metadata$CM7),1,0)

sfit <- survfit(Surv(tRFS, RFS)~CM1group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module1')

sfit <- survfit(Surv(tRFS, RFS)~CM2group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module2')

sfit <- survfit(Surv(tRFS, RFS)~CM3group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module3')

sfit <- survfit(Surv(tRFS, RFS)~CM4group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module4')

sfit <- survfit(Surv(tRFS, RFS)~CM5group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module5')

sfit <- survfit(Surv(tRFS, RFS)~CM6group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module6')

sfit <- survfit(Surv(tRFS, RFS)~CM7group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module7')


survival_metadata$event <- ifelse(survival_metadata$event=='Survival',0,1)
sfit <- survfit(Surv(tOS, event)~CM2group, data=survival_metadata)
ggsurvplot(data = survival_metadata,fit =  sfit,pval = T,xlim=c(0,2300),risk.table = T,x.break.by = 500,palette = c('#2ca9e1','#c53d43'))+ggtitle('Module2')

###### heatmap ---------
aaa <- prop.table(table(Normal_TIME$orig.ident,Normal_TIME$module),margin = 1) %>% as.data.frame()
bbb <- reshape2::dcast(data = aaa,formula = Var1~Var2,fill = 'Freq')
bbb <- tibble::column_to_rownames(bbb,var='Var1')
bbb[,c(1:ncol(bbb))] <- as.numeric(unlist(bbb[,c(1:ncol(bbb))]))

df <- bbb
library(pheatmap)
df$group <- apply(df, 1, function(x) names(df)[which.max(x)])
df$max_value_adjusted <- apply(df[,1:7], 1, function(x) {
  group <- names(df)[which.max(x)]
  max_value <- max(x)
  cm_number <- as.numeric(sub("CM", "", group))*10
  return(max_value + cm_number)
})
df_heatmap <- df[order(df$max_value_adjusted,decreasing = T),]
df_heatmap <- t(df_heatmap) %>% as.data.frame()
pheatmap::pheatmap(df_heatmap[1:7,],cluster_rows = F,cluster_cols = F,show_colnames = F,
                   border_color = NA,clustering_method = "ward.D2",color = colorRampPalette(colors = c("white","firebrick3",'firebrick4'))(100))


df <- bbb
df$Group <- apply(df, 1, function(x) names(x)[which.max(x)])
sort_df <- function(df, group_col) {
  df_list <- split(df, df[[group_col]])
  sorted_list <- lapply(names(df_list), function(cm) {
    df_cm <- df_list[[cm]]
    df_cm[order(-df_cm[, cm]), ]
  })
  do.call(rbind, sorted_list)
}
sorted_df <- sort_df(df, "Group")
df_heatmap <- sorted_df
df_heatmap$Group <- NULL
df_heatmap <- t(df_heatmap) %>% as.data.frame()
pheatmap::pheatmap(df_heatmap[1:7,],cluster_rows = F,cluster_cols = F,show_colnames = F,
                   border_color = NA,clustering_method = "ward.D2",color = colorRampPalette(colors = c("white","firebrick3",'firebrick4'))(100))

colors <- colorRampPalette(c("white", "firebrick3", "firebrick4"))(100)
breaks <- seq(min(df_heatmap), 1, length.out = length(colors) + 1)
pheatmap::pheatmap(df_heatmap[1:7, ], 
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_colnames = FALSE,
                   border_color = NA,
                   clustering_method = "ward.D2",
                   color = colors,
                   breaks = breaks)


