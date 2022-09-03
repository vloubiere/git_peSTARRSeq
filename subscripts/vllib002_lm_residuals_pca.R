setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh."]
pred <- readRDS("Rdata/CV_linear_model_vllib002.rds")$pred
dat[pred, diff:= log2FoldChange-i.predicted, on= c("L", "R")]

# Matrix for clustering
mat <- as.matrix(dcast(dat, L~R, value.var= "diff"), 1)
while(sum(is.na(mat))>0.025*length(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[, -which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}

# PCA
matL <- mat
matL <- apply(matL, 2, function(x) ifelse(is.na(x), mean(x, na.rm= T), x))
pcaL <- as.data.table(prcomp(matL)$x[, "PC1", drop= F], keep.rownames = "name")
matR <- t(mat)
matR <- apply(matR, 2, function(x) ifelse(is.na(x), mean(x, na.rm= T), x))
pcaR <- as.data.table(prcomp(matR)$x[, "PC1", drop= F], keep.rownames = "name")

# Order matri and create heatmap
x <- mat[order(-pcaL$PC1),]
pcaL <- pcaL[order(-PC1)]
x <- x[, order(pcaR$PC1)]
pcaR <- pcaR[order(PC1)]
cl <- vl_heatmap(x,
                 row_clusters = as.numeric(pcaL$PC1>0)+1,
                 col_clusters = as.numeric(pcaR$PC1>0)+1,
                 cluster_rows = F, 
                 cluster_cols = F,
                 breaks = c(-2,-0.25,0.25,2), 
                 col= c("cornflowerblue", "white", "white", "tomato"), 
                 legend_title = "Obs./Exp. (log2)", 
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
cl$rows[pcaL, PC1:= PC1, on= "name"]
cl$cols[pcaR, PC1:= PC1, on= "name"]
saveRDS(cl, 
        "Rdata/vllib002_lm_residuals_PCA.rds")

plot(cl)
