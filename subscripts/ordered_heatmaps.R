setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class== "enh./enh."]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
# Residuals
dat[, diff:= log2FoldChange-predict(model, new= dat)]
dat[, diff_add:= log2FoldChange-additive]
dat[, diff_mult:= log2FoldChange-multiplicative]
# Ordering
mean_median <- merge(unique(dat[, .(ID= L, median_L)]), 
                     unique(dat[, .(ID= R, median_R)]))
mean_median <- mean_median[, .(ID, rowMeans(.SD, na.rm= T)), .SDcols= c("median_L", "median_R")]
dat[mean_median, mean_median_L:= i.V2, on= "L==ID"]
dat[mean_median, mean_median_R:= i.V2, on= "R==ID"]

pca <- readRDS("Rdata/vllib002_PCA_LR.rds")
dat[as.data.table(pca$pca.L$x, keep.rownames = "ID"), PCA_L:= i.PC1, on= "L==ID"]
dat[as.data.table(pca$pca.R$x, keep.rownames = "ID"), PCA_R:= i.PC1, on= "R==ID"]

marg <- merge(unique(dat[, .(marg_L= mean(diff)), .(ID= L)]), 
              unique(dat[, .(marg_R= mean(diff)), .(ID= R)]))
marg[, mean_marg:= rowMeans(.SD, na.rm= T), .SDcols= c("marg_L", "marg_R")]
dat[marg, marg_L:= i.marg_L, on= "L==ID"]
dat[marg, marg_R:= i.marg_R, on= "R==ID"]
dat[marg, mean_marg_L:= mean_marg, on= "L==ID"]
dat[marg, mean_marg_R:= mean_marg, on= "R==ID"]

#-----------------------------------------------#
# Generate Heatmaps with different ordering ...
#-----------------------------------------------#
vars <- c("diff", "diff_add", "diff_mult")
ords <- list(c("median_L", "median_R"),
             c("mean_median_L", "mean_median_R"),
             c("PCA_L", "PCA_R"),
             c("marg_L", "marg_R"),
             c("mean_marg_L", "mean_marg_R"))

pdf("pdf/draft/ordered_residuals_heatmap.pdf", 
    width= 20, 
    height = 30)
par(mar= c(3.5,3.5,3,10),
    mgp= c(3,0.15,0),
    tcl= -0.2,
    mfrow= c(5,3))
for(ord in ords)
{
  for(var in vars)
  {
    # Matrix for clustering
    .f <- as.formula(paste0(ord[1], "+L~", ord[2], "+R"))
    mat <- as.matrix(dcast(dat, .f, value.var= var, sep = "__")[, -1], 1)
    while(sum(is.na(mat))>0.025*length(mat))
    {
      mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
      mat <- mat[, -which.max(apply(mat, 2, function(x) sum(is.na(x))))]
    }
    colnames(mat) <- unlist(tstrsplit(colnames(mat), "__", keep= 2))
    mat <- mat[rownames(mat) %in% colnames(mat),
               colnames(mat) %in% rownames(mat)]
    mat <- mat[nrow(mat):1,]
    
    # Plot
    if(grepl("mean", ord[1]))
      main <- paste(var, "residuals,", gsub("_L$", "", ord[1]), "ordering") else 
        main <- paste(var, "residuals,", ord[1], "*", ord[2], "ordering")
    cl <- vl_heatmap(mat,
                     breaks = c(-2,-0.25,0.25,2), 
                     cluster_rows= F,
                     cluster_cols= F,
                     col= c("cornflowerblue", "white", "white", "tomato"), 
                     legend_title = var, 
                     show_rownames = F,
                     show_colnames = F,
                     auto_margins = F,
                     main= main)
    ind_L <- dat[cl$rows[(order), name], median_L[1], L, on= "L"][match(rownames(mat), L), V1]
    ind_R <- dat[cl$cols[(order), name], median_R[1], R, on= "R"][match(colnames(mat), R), V1]
    
    # Left individual activities
    right <- par("usr")[1]-strwidth("M", cex= 0.5)
    width <- strwidth("M", cex= 3.5)
    rect(right-(ind_L/max(ind_L)*width), 
         rev(seq(ind_L))+0.5, 
         right, 
         rev(seq(ind_L))-0.5, 
         border= NA, 
         xpd= T, 
         col= adjustcolor("#0C3A0E", 0.7))
    ticks <- axisTicks(c(0, max(ind_L)), log= F)
    at <- right-(ticks/max(ind_L)*width)
    axis(3, 
         at = rev(range(at)), 
         labels = range(ticks),
         xpd= T, 
         line= 0.25, 
         cex.axis= 0.5)
    text(mean(at),
         par("usr")[4],
         "5' Individual\nact. (log2)",
         xpd= T, 
         pos= 3, 
         offset= 1.25,
         cex= 0.7)
    text(grconvertX(0.5, "line", "user"),
         mean(par("usr")[c(3,4)]),
         "5' enhancer",
         srt= 90,
         xpd= T)
    
    # Right individual activities
    par(mgp= c(3,0.35,0))
    top <- par("usr")[3]-strwidth("M", cex= 0.5)
    height <- -strheight("M", cex= 5)
    rect(seq(ind_R)+0.5,
         top, 
         seq(ind_R)-0.5,
         top+(ind_R/max(ind_R)*height),
         border= NA, 
         xpd= T, 
         col= adjustcolor("#0C3A0E", 0.7))
    ticks <- axisTicks(c(0, max(ind_R)), log= F)
    at <- top+(ticks/max(ind_R)*height)
    axis(2, 
         at = rev(range(at)), 
         labels = range(ticks),
         xpd= T, 
         line= 0.25, 
         cex.axis= 0.5,
         las= 2)
    text(par("usr")[1],
         mean(at),
         "3' Individual\nact. (log2)",
         xpd= T,
         pos= 2, 
         offset= 0.5,
         cex= 0.7)
    text(mean(par("usr")[c(1,2)]),
         grconvertY(0.5, "line", "user"),
         "3' enhancer",
         xpd= T)
  }
}
dev.off()






