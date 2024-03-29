setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")

# Matrix Inactive candidates
mat <- dcast(dat,
             -indL+L~indR+R,
             value.var = "residuals",
             sep = "__")
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
                 breaks = c(-2, 0, 2),
                 show_rownames= F,
                 show_colnames= F, 
                 legend_title= "Residuals (log2)",
                 plot= F,
                 legend.cex = 0.7)

# Plot
actLim <- 1
barLims <- c(-1,6)
barAxLims <- c(1,6)

pdf("pdf/draft/heatmap_residuals_ordered_ind_act.pdf",
    width = 3.85,
    height = 3.25)
layout(matrix(c(2,4,1,3), nrow= 2),
       widths= c(0.1,1),
       heights= c(1,0.1))
par(mar= c(0.2,0.2,0.2, 0.2),
    oma= c(3,3,2,6),
    mgp= c(2, 0.15, 0),
    tcl= -0.1,
    cex.axis= 0.6,
    las= 1,
    xpd= NA)
plot(hm)
left <- dat[rev(hm$rows$name), indL, on= "L", mult = "first"]
right <- dat[hm$cols$name, indR, on= "R", mult = "first"]
abline(h= min(which(left>actLim)),
       xpd= F)
abline(v= min(which(right>actLim)),
       xpd= F)
title(xlab= "3' activity",
      ylab= "5' activity")
par(lwd= 0.25)
legend(par("usr")[2]+strwidth("M")*0.75,
       par("usr")[4]-strheight("M")*9,
       c(NA, NA),
       fill= "lightgrey",
       border= T,
       xpd= NA,
       bty= "n",
       cex= 0.8)
legend(par("usr")[2]+strwidth("M")*0.75,
       par("usr")[4]-strheight("M")*9,
       c("Active", "Inactive"),
       dens= c(50, 0),
       border= F,
       xpd= NA,
       bty= "n",
       cex= 0.8)
# Left
par(yaxs= "i")
bar <- barplot(left,
               col= "lightgrey",
               space= 0,
               border= NA,
               horiz = T,
               xlim= barLims,
               xaxt= "n")
axis(side = 3,
     at = barAxLims,
     labels = barAxLims, 
     line = 0.2)
x <- left[left>actLim]
y <- bar[left>actLim]
polygon(c(0, x, 0),
        c(y[1], y, y[length(y)]),
        dens= 50)
# Right
par(yaxs= "r",
    xaxs= "i")
bar <- barplot(right,
               col= "lightgrey",
               space= 0,
               border= NA,
               ylim= barLims,
               yaxt= "n")
axis(side = 4,
     at = barAxLims,
     labels = barAxLims, 
     line = 0.2)
x <- bar[right>actLim]
y <- right[right>actLim]
polygon(c(x[1], x, x[length(x)]),
        c(0, y, 0),
        dens= 50)
dev.off()
