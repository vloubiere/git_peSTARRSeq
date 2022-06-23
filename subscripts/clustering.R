setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import dataset
if(!exists("feat"))
  feat <- readRDS("Rdata/final_300bp_enhancer_features_w_motifs.rds")
if(!exists("vl_screen"))
  vl_screen <- readRDS("Rdata/final_results_table.rds")

#-----------------------------------------------#
# Train activity based models and compute residuals
#-----------------------------------------------#
clean <- vl_screen[vllib=="vllib002" & class== "enh./enh."]
clean <- clean[L %in% feat[group %in% c("hk", "dev", "shared"), ID] 
               & R %in% feat[group %in% c("hk", "dev", "shared"), ID]]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
clean[, diff:= log2FoldChange-predict(model, new= clean)]

#-----------------------------------------------#
# Format data
#-----------------------------------------------#
dat <- merge(clean[, .(marg_L= mean(diff)), .(ID= L)],
             clean[, .(marg_R= mean(diff)), .(ID= R)])
dat[, PCC:= {
  .c <- merge(clean[.BY, .(ID= R, x= diff), on= "L==ID"],
              clean[.BY, .(ID= L, y= diff), on= "R==ID"])
  cor.test(.c$x, .c$y, method= "spearman")$estimate
}, ID]
plot(sort(dat$PCC))
setorderv(dat, "PCC", -1)
plot(merge(clean[dat[1, ID], .(ID= R, x= diff), on= "L"],
           clean[dat[1, ID], .(ID= L, y= diff), on= "R"])[, .(x, y)])
setorderv(dat, "PCC", 1)
plot(merge(clean[dat[1, ID], .(ID= R, x= diff), on= "L"],
           clean[dat[1, ID], .(ID= L, y= diff), on= "R"])[, .(x, y)])


#-----------------------------------------------#
# heatmap
#-----------------------------------------------#
mat <- as.matrix(dcast(clean, L~R, value.var= "diff"), 1)
while(sum(is.na(mat))>0.05*length(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[, -which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}

vl_heatmap(mat, 
           breaks = c(-2,-0.25,0.25,2), 
           cutree_rows = 4,
           cutree_cols = 4,
           col= c("cornflowerblue", "white", "white", "tomato"), 
           clustering_method = "ward.D")

# sym <- mat[, colnames(mat) %in% rownames(mat)]
# sym <- sym[rownames(mat) %in% colnames(mat),]
# apply(sym, )
# sym <- sym[order(rownames(sym)), order(colnames(sym))]
# vl_heatmap(mat, 
#            cluster_rows= F, 
#            cluster_cols= F, 
#            breaks = c(-2,-0.25,0.25,2),
#            col= c("cornflowerblue", "white", "white", "tomato"))