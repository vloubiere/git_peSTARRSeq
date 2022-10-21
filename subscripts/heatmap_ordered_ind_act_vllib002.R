setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")

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
Left <- unique(dat[rownames(mat), .(L, ctlL, actClassL, medianResidualsL), on= "L", keyby= indL])
Left[, smooth:= smooth.spline(medianResidualsL)$y]
Left[, col:= fcase(actClassL=="inactive", "lightgrey", 
                   default= "tomato"), actClassL]
Right <- unique(dat[colnames(mat), .(R, ctlR, actClassR, medianResidualsR), on= "R", keyby= indR])
Right[, smooth:= smooth.spline(medianResidualsR)$y]
Right[, col:= fcase(actClassR=="inactive", "lightgrey", 
                    default= "tomato"), actClassR]

# Sub matrix containing robust, usable active pairs
sub <- dat[L %in% Left[actClassL=="active" & smooth>0, L]
           & R %in% Right[actClassR=="active" & smooth>0, R]]
saveRDS(sub, "Rdata/vllib002_subset_active_pairs.rds")

# Plot
pdf("pdf/draft/heatmap_ordered_ind_act.pdf")
layout(matrix(c(7,2,8,5,1,4,6,3,9), nrow= 3),
       widths = c(1,10,1),
       heights = c(1,10,1))
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
abline(h= Left[actClassL=="inactive" | smooth>0, .N]+0.5)
abline(h= Left[col=="lightgrey", .N, col]$N+0.5)
abline(v= Right[actClassR=="inactive" | smooth>0, .N]+0.5)
abline(v= Right[col=="lightgrey", .N, col]$N+0.5)
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
          col= col,
          space= 0, 
          border= NA,
          horiz= T,
          xaxt= "n")
  axis(1,
       at= c(0,7),
       labels= c(0,7))
  par(mgp= c(2,0.4,0.2),
      mar= c(1,0.2,0.2,0.5))
  plot(medianResidualsL,
       y= .I,
       pch= 19,
       cex= 0.3,
       col= "grey80",
       xpd= NA,
       xaxt= "n",
       yaxt= "n",
       frame= F,
       xlab= NA,
       ylab= NA)
  lines(x= smooth.spline(medianResidualsL)$y,
        y= .I,
        col= "grey80",
        lwd= 2)
  axis(1,
       at= c(-1, 0.5),
       labels= c(-1, 0.5))
  abline(v= 0, lty= 2)
}]
# Right individual act
Right[, {
  par(mgp= c(2,0.5,0.2),
      mar= c(0.5,1,0.2,0.2))
  barplot(indR,
          col= col,
          space= 0, 
          border= NA,
          yaxt= "n")
  axis(2,
       at= c(0, 7),
       labels= c(0, 7))
  par(mgp= c(2,0.5,0.2),
      mar= c(0.5,1,0.5,0.2))
  plot(medianResidualsR,
       pch= 19,
       cex= 0.3,
       col= "grey80",
       xpd= NA,
       xaxt= "n",
       yaxt= "n",
       frame= F,
       xlab= NA,
       ylab= NA)
  lines(smooth.spline(medianResidualsR)$y,
        col= "grey80",
        lwd= 2)
  axis(2,
       at= c(-0.8, 0.5),
       labels= c(-0.8, 0.5))
  abline(h= 0, lty= 2)
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



