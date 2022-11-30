setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_dev_pairs_vllib002_with_predictions.rds")

# Ordering
dat[, L:= factor(L, unique(L[order(indL, decreasing = T)]))]
dat[, R:= factor(R, unique(R[order(indR)]))]

# Matrix
Cc <- circlize::colorRamp2(c(-2,0,2), col= c("blue", "grey", "red"))
mat <- dcast(dat, L~R, value.var = "residuals")
mat <- as.matrix(mat, 1)
while(sum(is.na(mat))>0.1*nrow(mat)*ncol(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}

# Left and right margins object
Left <- unique(dat[rownames(mat), .(L, ctlL, actClassL, meanResidualsL), on= "L", keyby= indL])
Left[, smooth:= smooth.spline(meanResidualsL)$y]
Right <- unique(dat[colnames(mat), .(R, ctlR, actClassR, meanResidualsR), on= "R", keyby= indR])
Right[, smooth:= smooth.spline(meanResidualsR)$y]

# Plot
pdf("pdf/draft/heatmap_ordered_ind_act.pdf")
layout(matrix(c(7,2,8,5,1,4,6,3,9), nrow= 3),
       widths = c(0.8,10,1.2),
       heights = c(1.2,10,0.8))
par(mar= c(1,1,0.2,0.2),
    tcl= -0.2,
    las= 1,
    col.axis= "black",
    xaxs= "i",
    yaxs= "i")
# Heatmap
hm <- vl_heatmap(mat,
                 cluster_rows= F,
                 cluster_cols= F,
                 breaks = c(-1,-0.25, 0.25, 1), 
                 col= c("cornflowerblue", "white", "white", "tomato"),
                 show_rownames = F,
                 show_colnames = F, 
                 show_legend = F)
# coordinates
left <- grconvertX(0, "nfc", "user")
lw <- diff(grconvertX(c(0,1), "line", "user"))
bottom <- grconvertY(0, "nfc", "user")
lh <- diff(grconvertY(c(0,1), "line", "user"))
# Left individual act
lim <- c(-1.5, 8.5)
Left[, {
  par(mgp= c(2,0.4,0.2),
      mar= c(1,0.5,0.2,0.2))
  barplot(indL,
          space= 0, 
          border= NA,
          horiz= T,
          xaxt= "n")
  axis(1,
       at= c(0,4),
       labels= c(0,4))
  par(mgp= c(2,0.4,0.2),
      mar= c(1,0.2,0.2,0.5))
  barplot(meanResidualsL,
          space= 0, 
          border= NA,
          horiz= T,
          xaxt= "n")
  abline(v= 0, lwd= 0.5)
  axis(1,
       at= c(-1, 0, 1),
       labels= c(-1, 0, 1))
}]
# Right individual act
Right[, {
  par(mgp= c(2,0.5,0.2),
      mar= c(0.5,1,0.2,0.2))
  barplot(indR,
          space= 0, 
          border= NA,
          yaxt= "n")
  axis(2,
       at= c(0, 4),
       labels= c(0, 4))
  par(mgp= c(2,0.5,0.2),
      mar= c(0.5,1,0.5,0.2))
  barplot(meanResidualsR,
          space= 0, 
          border= NA,
          yaxt= "n")
  axis(2,
       at= c(-1, 0, 0.6),
       labels= c(-1, 0, 0.6))
  abline(h= 0, lwd= 0.5)
}]
par(mar= c(0,0,0,0))
plot.new()
vl_heatkey(breaks = hm$breaks, 
           col= hm$col, 
           left = 0.1, 
           top = 0.7, 
           height = 0.6, 
           width = 0.1, 
           main = "Obs./Exp. (log2)", 
           main.cex = 0.6)
dev.off()