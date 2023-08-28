setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat <- dat[actL=="Inactive" & actR=="Inactive"]

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
actLim <- 0
barLims <- c(1,-1)
barAxLims <- c(0,1)

left <- dat[rev(hm$rows$name), indL, on= "L", mult = "first"]
right <- dat[hm$cols$name, indR, on= "R", mult = "first"]
Cc <- c("grey60", "#74C27A", "grey30", "grey30", "tomato", "grey30")
leftType <- factor(unlist(tstrsplit(rev(hm$rows$name), "_", keep= 1)))
rightType <- factor(unlist(tstrsplit(hm$cols$name, "_", keep= 1)))

pdf("pdf/draft/heatmap_residuals_ordered_ind_act_inactive_pairs.pdf",
    width = 3.4,
    height = 2.7)
layout(matrix(c(c(3,6,7,2,8,9,1,4,5)), nrow= 3),
       widths= c(0.15,0.06,1),
       heights= c(1,0.06,0.15))
par(mar= c(0.2,0.2,0.2,0.2),
    oma= c(3,3,1,7),
    mgp= c(2, 0.15, 0),
    tcl= -0.1,
    cex.axis= 0.6,
    las= 1,
    xpd= NA)
plot(hm)
abline(h= min(which(left>actLim)),
       xpd= F)
abline(v= min(which(right>actLim)),
       xpd= F)
title(xlab= "3' activity",
      ylab= "5' activity",
      line = 3)
par(lwd= 0.25)
legend(par("usr")[2]+strwidth("M")*0.75,
       par("usr")[4]-strheight("M")*9,
       legend = c("Random seq.", "Inactive", "Developmental", "Housekeeping"),
       fill= c("grey60", "grey30", "#74C27A", "tomato"),
       border= T,
       xpd= NA,
       bty= "n",
       cex= 0.8)
legend(par("usr")[2]+strwidth("M")*0.75,
       par("usr")[4]-strheight("M")*15,
       c(NA, NA),
       fill= "lightgrey",
       border= T,
       xpd= NA,
       bty= "n",
       cex= 0.8)
legend(par("usr")[2]+strwidth("M")*0.75,
       par("usr")[4]-strheight("M")*15,
       c("Active", "Inactive"),
       dens= c(50, 0),
       border= F,
       xpd= NA,
       bty= "n",
       cex= 0.8)
# Left
# Left
par(yaxs= "i")
plot.new()
rasterImage(matrix(Cc[leftType], ncol= 1),
            xleft = 0,
            ybottom = 1,
            xright = 1,
            ytop = 0,
            interpolate = F)
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
plot.new()
rasterImage(matrix(Cc[rightType], nrow= 1),
            xleft = 0,
            ybottom = 0,
            xright = 1,
            ytop = 1,
            interpolate = F)
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