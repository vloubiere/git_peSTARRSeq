setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
dat[, L:= factor(L, unique(L[order(indL, decreasing = T)]))]
dat[, R:= factor(R, unique(R[order(indR)]))]

Cc <- circlize::colorRamp2(c(-2,0,2), col= c("blue", "grey", "red"))
mat <- dcast(dat, L~R, value.var = "residuals")
mat <- as.matrix(mat, 1)
while(sum(is.na(mat))>0.05*nrow(mat)*ncol(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}

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
                 breaks = c(-1.25,-0.25, 0.25, 1.25), 
                 col= c("cornflowerblue", "white", "white", "tomato"),
                 show_rownames = F,
                 show_colnames = F, 
                 show_legend = F)
# Margins
lim <- c(-1.5, 8.5)
norm.score <- function(x, all)
  (x-min(all))/diff(range(all))*diff(lim)+min(lim)
# Left and right objects
Left <- unique(dat[rownames(mat), .(L, ctlL, actClassL, meanResidualsL), on= "L", keyby= indL])
Left[, col:= fcase(ctlL, "grey20",
                   actClassL=="active", "tomato",
                   default= "grey80")]
Right <- unique(dat[colnames(mat), .(R, ctlR, actClassR, meanResidualsR), on= "R", keyby= indR])
Right[, col:= fcase(ctlR, "grey20",
                    actClassR=="active", "tomato",
                    default= "grey80")]
left <- grconvertX(0, "nfc", "user")
lw <- diff(grconvertX(c(0,1), "line", "user"))
bottom <- grconvertY(0, "nfc", "user")
lh <- diff(grconvertY(c(0,1), "line", "user"))
Left[, .N, rleid(actClassL)][, {
  rect(left+0.1*lw, 
       c(1, cumsum(N)[-2])-0.5, 
       left+0.9*lw, 
       cumsum(N)+0.5, 
       border= NA, 
       xpd= T,
       col= c("lightgrey", "tomato"))
  text(left+0.5*lw,
       cumsum(N)-N/2,
       c("Inactive", "Active"),
       srt= -90,
       xpd= NA)
  abline(h= N[1]+0.5)
}]
Right[, .N, rleid(actClassR)][, {
  rect(c(1, cumsum(N)[-2])-0.5, 
       bottom+lh*0.1, 
       cumsum(N)+0.5, 
       bottom+lh*0.9, 
       border= NA, 
       xpd= T,
       col= c("lightgrey", "tomato"))
  text(cumsum(N)-N/2,
       bottom+lh*0.5,
       c("Inactive", "Active"),
       xpd= NA)
  abline(v= N[1]+0.5)
}]
# Left Frame
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
  plot(meanResidualsL,
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
  lines(x= smooth.spline(meanResidualsL)$y,
        y= .I,
        col= "grey80",
        lwd= 2)
  axis(1,
       at= c(-1.5, 0.5),
       labels= c(-1.5, 0.5))
  abline(v= 0, lty= 2)
}]
# Right Frame
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
  plot(meanResidualsR,
       pch= 19,
       cex= 0.3,
       col= "grey80",
       xpd= NA,
       xaxt= "n",
       yaxt= "n",
       frame= F,
       xlab= NA,
       ylab= NA)
  lines(smooth.spline(meanResidualsR)$y,
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