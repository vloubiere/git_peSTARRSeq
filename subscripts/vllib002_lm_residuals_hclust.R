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
cl <- vl_heatmap(mat,
                 breaks = c(-2,-0.25,0.25,2), 
                 cutree_rows = 4,
                 cutree_cols = 4,
                 col= c("cornflowerblue", "white", "white", "tomato"), 
                 # clustering_method = "ward.D2",
                 legend_title = "Obs./Exp. (log2)", 
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
cl$rows[dat, median:= median_L, on= "name==L", mult= "first"]
cl$cols[dat, median:= median_R, on= "name==R", mult= "first"]
saveRDS(cl, "Rdata/vllib002_lm_residuals_hclust.rds")

plot(cl)
