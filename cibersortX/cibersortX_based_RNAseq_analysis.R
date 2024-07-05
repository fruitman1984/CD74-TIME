library(tidyvers)
library(pheatmap)
library(survival)
library(survminer)

cibersort <- read.csv('CIBERSORTx_Job91_Results.csv',header = T,row.names = 1) %>% as.data.frame()

cibersort$CM1 <- rowSums(cibersort[,c('CD4T_CD52','Mac_CD74_MT1X','Plasma_IgM_IgA_IgG','CD8T_HAVCR2','CD8T_MKI67','CD8T_TYMS')])
cibersort$CM2 <- rowSums(cibersort[,c('Mac_TLR4','Mac_IL1B','cDC_CD1C','Plasma_IgA_IgG_SDC1','Plasma_IgA_IgG_CD27','Plasma_IgA_IgG_MKI67','Mono_IFI30','Mono_MPO')])
cibersort$CM3 <- rowSums(cibersort[,c('NK_CREM','B_IgM_IgD_FOXP1','Mono_CCNL1','NK_NKTR','CD4T_SF1','CD8T_FGFBP2','CD8T_CD226', 'NK_IKZF2','Plasma_IgA_IgG_IGLC1')])
cibersort$CM4 <- rowSums(cibersort[,c('CD8T_CXCR4','CD8T_CD160','NK_CD56','B_CD24_IER2','NK_CCL3','NK_IFNG', 'B_IgM_IgD_CD69','Mono_JUN','Mono_VCAN','Mono_PLBD1')])
cibersort$CM5 <- rowSums(cibersort[,c('Mac_LY6E', 'B_IgM_IgD_FCRL1','NK_IFI6','CD8T_IFIT3','B_CD34','pDC_LILRA4',
                                      'pDC_JCHAIN','B_IgM_IgD_CD24_CD1c', 'Plasma_IgM_IgA_CD27','NK_CD16','NK_HLADRA','CD4T_ANXA1','NK_PDCD4','CD8T_CCL5')])
cibersort$CM6 <- rowSums(cibersort[,c('CD4T_FOXP3','B_CD11c_HLADR','CD4T_LTB', 'NK_GZMH','gdT_TRGC1_TRDC','CD8T_CD16','Mono_S100A10','Mac_CD74_CD49d','MAIT_SLC4A10','B_CD24_TXNIP','CD8T_LYAR','B_CD11b_ANXA2')])
cibersort$CM7 <- rowSums(cibersort[,c('NK_MT1X','Mono_APOC1','Mac_CD163','Neu_CD66b','Mono_CYP11B1','Neu_CD16b',
                                      'CD4T_MAL','CD8T_LEF1','NK_PRF1','B_IgM_IgD_TCL1A','NK_ITGB1','NK_ADGRG1')])

cibersort2 <- cibersort[,paste0('CM',1:7)]

survival <- read.csv("E:/MyData/AMLmap/meta.csv",header = T,row.names = 1)
samples <- intersect(rownames(cibersort2),rownames(survival))
cibersort2 <- cibersort2[samples,]
survival <- survival[samples,]
survival_ciber <- cbind(cibersort2,survival)


survival_ciber$event <- ifelse(survival_ciber$event == 'Alive',0,1)
survival_ciber$OS <- as.numeric(round(survival_ciber$survival_years*365))

survival_ciber$cm1group <- ifelse(survival_ciber$CM1 > quantile(survival_ciber$CM1,0.6),1,0)
fit <- survfit(Surv(OS, event) ~ cm1group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3
) 
#### 5*7

survival_ciber$cm1group <- ifelse(survival_ciber$CM1 >= quantile(survival_ciber$CM1,0.498),1,0)
fit <- survfit(Surv(OS, event) ~ cm1group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3
)

survival_ciber$cm1group <- ifelse(survival_ciber$CM1 >= quantile(survival_ciber$CM1,0.6),1,0)
fit <- survfit(Surv(OS, event) ~ cm1group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3
) 

survival_ciber$cm2group <- ifelse(survival_ciber$CM2 > quantile(survival_ciber$CM2,0.51),1,0)
fit <- survfit(Surv(OS, event) ~ cm2group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3
) 

survival_ciber$cm2group <- ifelse(survival_ciber$CM2 > quantile(survival_ciber$CM2,0.458),1,0)
fit <- survfit(Surv(OS, event) ~ cm2group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3
) 

survival_ciber$cm3group <- ifelse(survival_ciber$CM3 > quantile(survival_ciber$CM3,0.49),1,0)
fit <- survfit(Surv(OS, event) ~ cm3group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm4group <- ifelse(survival_ciber$CM4 > quantile(survival_ciber$CM4,0.495),1,0)
fit <- survfit(Surv(OS, event) ~ cm4group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm4group <- ifelse(survival_ciber$CM4 > quantile(survival_ciber$CM4,0.7),1,0)
fit <- survfit(Surv(OS, event) ~ cm4group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm5group <- ifelse(survival_ciber$CM5 > quantile(survival_ciber$CM5,0.493),1,0)
fit <- survfit(Surv(OS, event) ~ cm5group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm5group <- ifelse(survival_ciber$CM5 > quantile(survival_ciber$CM5,0.7),1,0)
fit <- survfit(Surv(OS, event) ~ cm5group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm6group <- ifelse(survival_ciber$CM6 > quantile(survival_ciber$CM6,0.493),1,0)
fit <- survfit(Surv(OS, event) ~ cm6group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm6group <- ifelse(survival_ciber$CM6 > quantile(survival_ciber$CM6,0.493),1,0)
fit <- survfit(Surv(OS, event) ~ cm6group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm6group <- ifelse(survival_ciber$CM6 > quantile(survival_ciber$CM6,0.75),1,0)
fit <- survfit(Surv(OS, event) ~ cm6group, data = survival_ciber)

survival_ciber$cm7group <- ifelse(survival_ciber$CM7 > quantile(survival_ciber$CM7,0.501),1,0)
fit <- survfit(Surv(OS, event) ~ cm7group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

survival_ciber$cm7group <- ifelse(survival_ciber$CM7 > quantile(survival_ciber$CM7,0.63),1,0)
fit <- survfit(Surv(OS, event) ~ cm7group, data = survival_ciber)
ggsurvplot(fit, pval = TRUE,palette=c("#2ca9e1", "#c53d43"),title="Overall survival",
           risk.table = T, risk.table.height = 0.3)

####### heatmap --------
df <- cibersort2[,c('CM1','CM2','CM3','CM4','CM5','CM6','CM7')] %>% as.data.frame()
df$CM1 <- df$CM1 *1.6
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
df_heatmap$CM1 <- df_heatmap$CM1/1.6
df_heatmap <- t(df_heatmap) %>% as.data.frame()

anno_col <- read.csv("E:/MyData/AMLmap/meta.csv",header = T,row.names = 1)
anno_col <- anno_col[,'study',drop=F]

annotation_colors <- list(study=c('100LUMC_AML'='#8b968d',"BEAT_AML"='#c5c56a',"LEUCEGENE_AML"='#5c9291',"TARGET_AML"='#946243',"TCGA_AML"='#f7b977'))

pheatmap::pheatmap(df_heatmap,cluster_rows = F,cluster_cols = F,show_colnames = F,annotation_col = anno_col,annotation_colors = annotation_colors,
                   border_color = NA,color = colorRampPalette(colors = c("white","firebrick3",'firebrick4'))(101))