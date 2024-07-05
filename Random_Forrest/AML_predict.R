setwd('/GPUFS/sysu_zhaom_1/QG/300/predict')
load('/GPUFS/sysu_zhaom_1/QG/300/seurat_data_merge.RData')

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(randomForest)


## load classifier
load("/GPUFS/sysu_zhaom_1/QG/300/try/Train/RF1.RData", verbose = TRUE)   # as generated using the RFtrain script
load("/GPUFS/sysu_zhaom_1/QG/300/try/Train/RF2.RData", verbose = TRUE)

# X <- log(E+1)
X <- as.matrix(GetAssayData(seurat_data_merge,assay = 'RNA',slot = 'data'))
## predict, within seconds
# RF1

#X.predict.all <- predict(rf.all.inner, newdata = t(X[rownames(rf.all.inner$importance), ]), type="prob")

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
#X.predict.tumor <- predict(rf.tumor.inner, newdata = t(X[rownames(rf.tumor.inner$importance), ]), type="prob")

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