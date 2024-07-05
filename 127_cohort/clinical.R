######## clinical ------
load('/GPUFS/sysu_zhaom_1/QG/CD74TIME/RData/Normal_TIME.RData')
setwd('/GPUFS/sysu_zhaom_1/QG/CD74TIME/CD74TIME/module')
mutation <- read.csv('mutation_3group_127.csv',header = T)
survial_meta <- read.csv('survival_metadata.csv',header = T)
clinical <- read.csv('clinical.csv',header = T)

mutation <- subset(mutation,Patient %in% unique(survial_meta$X))
clinical <- subset(clinical,Patient %in% unique(survial_meta$X))
survial_meta <- tibble::column_to_rownames(survial_meta,var = 'X') %>% t() %>% as.data.frame()
clinical <- tibble::column_to_rownames(clinical,var = 'Patient') %>% t() %>% as.data.frame()
mutation_dcast <- reshape2::dcast(mutation,formula = Patient~Gene,value.var = 'Mutation', fun.aggregate = max,fill = 0)
mutation_dcast <- tibble::column_to_rownames(mutation_dcast,var = 'Patient')
for (i in 1:(ncol(mutation_dcast))) {
  mutation_dcast[,i] <- ifelse(mutation_dcast[,i] == 0, 0,1)
}

aaa <- prop.table(table(Normal_TIME$orig.ident,Normal_TIME$module),margin = 1) %>% as.data.frame()
bbb <- reshape2::dcast(data = aaa,formula = Var1~Var2,fill = 'Freq')
bbb <- tibble::column_to_rownames(bbb,var='Var1')
bbb[,c(1:ncol(bbb))] <- as.numeric(unlist(bbb[,c(1:ncol(bbb))]))
df <- bbb
df$CM1 <- df$CM1 * 1.5
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
df_heatmap$CM1 <- df_heatmap$CM1/1.5
df_heatmap <- t(df_heatmap) %>% as.data.frame()

clinical <- clinical[,colnames(df_heatmap)]
survial_meta <- survial_meta[,colnames(df_heatmap)]
mutation_dcast <- mutation_dcast[colnames(df_heatmap),]
rownames(mutation_dcast) <- colnames(df_heatmap)
mutation_dcast <- t(mutation_dcast) %>% as.data.frame()
mutation_dcast[is.na(mutation_dcast)] <- 3
mutation_dcast[,c(1:ncol(mutation_dcast))] <- as.numeric(unlist(mutation_dcast[,c(1:ncol(mutation_dcast))]))

library(ComplexHeatmap)

FAB_df <- clinical['FAB',]
FAB_df <- apply(FAB_df, c(1, 2), function(x) paste0("M", x))
col_FAB <- c(M0='#669966',M1='#66CCFF',M2='#FFCC00',M4='#9966FF',M5='#FF6633')
ht1 = Heatmap(FAB_df, name = "FAB", col = col_FAB, row_title = "FAB",
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht1

ELN2022_df <- clinical['ELN2022',]
col_ELN2022 <- c(Favorable='#F9E79F',Intermediate='#F39C12',Adverse='#D35400',Unknown='white')
ht2 = Heatmap(ELN2022_df, name = "ELN2022", row_title = "ELN2022", col = col_ELN2022,
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht2

Age_df <- clinical['Age',]
Age_df[,c(1:ncol(Age_df))] <- as.numeric(unlist(Age_df[,c(1:ncol(Age_df))]))
col_AGE <- colorRamp2(c(13,83), colors = c('#EAFAF1','#196F3D'))
ht3 = Heatmap(Age_df, name = "AGE", row_title = "AGE", col = col_AGE,
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht3

Gender_df <- clinical['Gender',]
Gender_df <- apply(Gender_df, c(1, 2), function(x) ifelse(x==1,'Male','Female'))
col_Gender <- c('Male'='#7499F4',Female='#F4749C')
ht4 = Heatmap(Gender_df, name = "Gender", row_title = "Gender", col = col_Gender,
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht4

WBC_df <- clinical['WBC',]
WBC_df[,c(1:ncol(WBC_df))] <- as.numeric(unlist(WBC_df[,c(1:ncol(WBC_df))]))
ht5 = Heatmap(WBC_df, name = "WBC", row_title = "WBC", col = c('#EAF2F8','#5499C7','#1F618D','#1F618D','#1F618D'),
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht5


RR_df <- survial_meta['RR',]
col_RR <- c(remission='#B6D4FA',relapsed='#ECB6FA',refractory='#FAEAB6',Death_without_remission='#FFABAB')
ht6 = Heatmap(RR_df, name = "First Event", row_title = "First Event", col = col_RR,
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht6

tREL <- survial_meta['tREL',]
tREL[,c(1:ncol(tREL))] <- as.numeric(unlist(tREL[,c(1:ncol(tREL))]))
ht7 = Heatmap(tREL, name = "RFS", row_title = "RFS", col =  c('#FFE2E2','#D10505'),
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht7

FG_df <- clinical[c("CBFBMYH11","AML1ETO"),]
FG_df <- apply(FG_df, c(1, 2), function(x) ifelse(x==1,'MUT',
                                                  ifelse(x==0,'WT',
                                                         ifelse(x==2,'Unknown','Others'))))
mutation_df <- mutation_dcast[c('TET2','FLT3','CEBPA','NPM1','ASXL1','DNMT3A','GATA2','TET1','CD101','FAT1','IDH2','DDX18','NRAS','KIT','RUNX1'),]
mutation_df <- apply(mutation_df, c(1, 2), function(x) ifelse(x==1,'MUT',
                                                              ifelse(x==0,'WT',
                                                                     ifelse(x==3,'Unknown','Others'))))
mut_df <- rbind(FG_df,mutation_df)
ht8 = Heatmap(mut_df, name = "Mutation", row_title = "Mutation", col =  c(WT='lightgrey',MUT='black',Unknown='white'),
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht8

colors <- colorRampPalette(c("white", "firebrick3", "firebrick4"))(100)
breaks <- seq(min(df_heatmap), 1, length.out = length(colors) )
colors <- colorRamp2(breaks, colors)
ht10086 = Heatmap(df_heatmap, name = "Module", col = colors, row_title = "Module",
                  rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht10086

ht_list = ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% ht7 %v% ht8 %v% ht10086
draw(ht_list)
### 10*7
matrix_clin <- rbind(FAB_df,ELN2022_df,Age_df,Gender_df,WBC_df,RR_df,tREL,mut_df,df_heatmap)


matrix_clin_t <- t(matrix_clin) %>% as.data.frame()
last_seven_colnames <- colnames(matrix_clin_t)[(ncol(matrix_clin_t)-6):ncol(matrix_clin_t)]
matrix_clin_t$group <- apply(matrix_clin_t[, last_seven_colnames], 1, function(x) names(x)[which.max(x)])
matrix_clin_t$group[1:8] <- 'CM1'
group_df <- as.data.frame(t(matrix_clin_t[,'group',drop=F]))
write.csv(matrix_clin_t,file = 'matrix_clin.csv')
ht0 = Heatmap(group_df, name = "Module", row_title = "Module", col =  c(CM1="#FB8072",CM2="#BEBADA",CM3="#FFFFB3",CM4="#8DD3C7",CM5="#FDB462",CM6="#B3DE69",CM7="#80B1D3"),
              rect_gp = gpar(col = 'white',lwd=2),cluster_columns = F,cluster_rows = F, show_column_names = FALSE)
ht0
ht_list = ht0 %v% ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% ht7 %v% ht8 %v% ht10086
draw(ht_list)

#### exam
matrix_exam <- rbind(FAB_df,ELN2022_df,Age_df,Gender_df,WBC_df,RR_df,tREL,mut_df) %>% t() %>% as.data.frame()
matrix_exam$group <- matrix_clin_t$group

library(dplyr)
library(broom)

matrix_exam$Age <- as.numeric(matrix_exam$Age)
matrix_exam$WBC <- as.numeric(matrix_exam$WBC)
matrix_exam$tREL <- as.numeric(matrix_exam$tREL)

test_column <- function(data, column) {
  data_subset <- data[, c(column, "group"), drop = FALSE]
  data_subset <- data_subset[data_subset[[column]] != 'Unknown', ]
  data_subset[[column]] <- factor(data_subset[[column]])
  if(is.numeric(data_subset[[column]])) {
    anova_result <- aov(data_subset[[column]] ~ group, data = data_subset)
    return(tidy(anova_result))
  } else if(is.factor(data_subset[[column]])) {
    contingency_table <- table(data_subset[[column]], data_subset$group)
    if(all(chisq.test(contingency_table)$expected >= 5)) {
      chisq_result <- chisq.test(contingency_table)
      return(tidy(chisq_result))
    } else {
      fisher_result <- fisher.test(contingency_table, workspace = 2e8)
      return(tidy(fisher_result))
    }
  } else {
    return(data.frame())
  }
}

results <- lapply(names(df)[-which(names(df) == "group")], function(column) {
  test_column(df, column)
})

results