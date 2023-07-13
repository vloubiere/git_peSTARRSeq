setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")

# Make matrix and clean
mat <- dcast(dat, -indL+L~indR+R, value.var = "log2FoldChange", sep = "__")
mat <- as.matrix(mat[, -1], 1)
colnames(mat) <- unlist(tstrsplit(colnames(mat), "__", keep= 2))
while(sum(is.na(mat))>0.05*nrow(mat)*ncol(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}
hm <- vl_heatmap(mat, 
                 cluster_rows= F,
                 cluster_cols= F, 
                 breaks = seq(-2, 8, 2),
                 col= viridis::viridis(6), 
                 show_rownames= F,
                 show_colnames= F, 
                 legend_title= "Activity (log2)")

Cc <- circlize::colorRamp2(hm$breaks, hm$col)
barplot_width_cex <- 1.5

pdf("pdf/draft/heatmap_ordered_ind_act.pdf", width = 5, height = 5)
par(mar= c(5.1, 4.1, 4.1, 5.5),
    mgp= c(2, 0.15, 0),
    tcl= -0.1,
    cex.axis= 0.6,
    las= 1)
plot(hm)
# Add left and right barplots
for(side in c("left", "right"))
{
  var <- if(side=="left")
    rev(dat[L %in% rownames(mat)][order(-indL), unique(indL)]) else
      dat[R %in% colnames(mat)][order(indR), unique(indR)]
  width <- if(side=="left") # Line width
    diff(grconvertX(c(0, 1), "lines", "user")) else
      diff(grconvertY(c(0, 1), "lines", "user"))
  minPos <- par("usr")[ifelse(side=="left", 1, 3)]-width*0.5 # min value position
  .center <- minPos-(0-min(var))/diff(range(var))*width*barplot_width_cex
  # Add axis
  .leg <- c(1, floor(max(var)))
  axis(ifelse(side=="left", 3, 4),
       at= minPos-(.leg-min(var))/diff(range(var))*width*barplot_width_cex, 
       labels = .leg,
       xpd= T, 
       line = 0.25)
  # Compute barplots
  coor <- data.table(x= minPos-(var-min(var))/diff(range(var))*width*barplot_width_cex,
                     y= seq(var),
                     act= var>=1)
  # Full polygon
  full <- coor[, .(x= c(.center, x, .center),
                   y= c(y[1], y, y[.N]))]
  if(side!="left")
    setnames(full, c("x", "y"), c("y", "x"))
  full[, polygon(x, 
                 y,
                 col= "lightgrey",
                 border= NA,
                 lwd= 0.25,
                 xpd= NA)]
  # Active enhancers
  act <- coor[(act), .(x= c(.center, x, .center),
                       y= c(y[1], y, y[.N]))]
  if(side!="left")
    setnames(act, c("x", "y"), c("y", "x"))
  act[, {
    polygon(x, 
            y, 
            dens= 50,
            lwd= 0.25,
            xpd= NA)
  }]
  # Add line
  if(side=="left")
    abline(h= min(which(var>=1)), col= "white") else
      abline(v= min(which(var>=1)), col= "white")
}
# Global labels
title(xlab= "3' activity",
      ylab= "5' activity")
par(lwd= 0.25)
legend(par("usr")[2],
       mean(par("usr")[c(3,4)]),
       c(NA, NA),
       fill= "lightgrey",
       border= T,
       xpd= NA,
       bty= "n",
       cex= 0.8)
legend(par("usr")[2],
       mean(par("usr")[c(3,4)]),
       c("Active", "Inactive"),
       dens= c(50, 0),
       border= F,
       xpd= NA,
       bty= "n",
       cex= 0.8)
# Add left and right act heatmaps
# left <- matrix(dat[L %in% rownames(mat)][order(-indL), unique(indL)], ncol= 1)
# rasterImage(Cc(left), 
#             par("usr")[1]-strwidth("M")*2.5,
#             par("usr")[3],
#             par("usr")[1]-strwidth("M"),
#             par("usr")[4], 
#             xpd= NA)
# Add right act heatmap
# right <- matrix(dat[R %in% colnames(mat)][order(indR), unique(indR)], nrow= 1)
# rasterImage(Cc(right), 
#             par("usr")[1],
#             par("usr")[3]-strheight("M")*2.5,
#             par("usr")[2],
#             par("usr")[3]-strheight("M"), 
#             xpd= NA)
dev.off()
