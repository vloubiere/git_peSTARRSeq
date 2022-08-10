setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
options(scipen= 1)
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
dat <- model$vllib002_dat
# Select examples
ex <- dat[list("dev_strong_B_00302",
               c("dev_medium_C_00480", "dev_medium_C_00466", "dev_medium_B_00585")), on= c("L", "R")]
ex[, col:= c("tomato", "limegreen", "cornflowerblue")]

#-----------------------------------------------#
# Plot
#-----------------------------------------------#
Cc <- colorRampPalette(c("white", "grey90", "grey60", "grey20"))
pdf("pdf/draft/Figure_2A.pdf",
    height = 6,
    width = 6)
layout(matrix(c(1,2,2,2,
                3,4,4,11,
                3,5,8,11,
                3,6,9,11,
                3,7,10,11,
                3,12,13,14), 
              ncol= 4, 
              byrow = T),
       widths = c(1,0.2,0.2,0.6),
       heights =  c(1, 0.122, rep(0.256, 3), 0.11))
# Smoothscatter add, mult and lm
par(las= 1,
    mar= c(3,4,2,0.5),
    mgp= c(1.75, 0.5, 0),
    tcl= -0.2,
    cex= 1)
smoothScatter(dat[, .(`Expected additive (log2)`= additive, 
                      `Activity (log2)`= log2FoldChange)],
              main= "Additive",
              colramp = Cc)
abline(0,1,lty=2)
smoothScatter(dat[, .(`Expected multiplicative (log2)`= multiplicative, 
                      `Activity (log2)`=log2FoldChange)],
              main= "Multiplicative",
              colramp = Cc)
abline(0,1,lty=2)
smoothScatter(dat[, .(`Predicted (log2)`= predicted, 
                      `Activity (log2)`= log2FoldChange)],
              main= "Linear model",
              colramp = Cc)
abline(0,1,lty=2)
ex[, points(predicted, 
            log2FoldChange, 
            col= col, 
            pch= 16), .(L, R, col)]
legend(par("usr")[1]-strwidth("M", cex= 0.6),
       par("usr")[4],
       c(model$eq, 
         paste0("R²= ", round(model$rsq$Rsquare, 2))),
       bty= "n",
       cex= 0.6)

# Examples barplot
par(mar= c(0,0,0,0),
    lwd= 0.1)
plot.new()
segments(0.3, 
         seq(0.3, 0.7, length.out=3), 
         0.9,
         seq(0.3, 0.7, length.out=3))
rect(rep(c(0.375, 0.675), each= 3),
     rep(seq(0.3, 0.7, length.out=3)-0.06, 2),
     rep(c(0.375, 0.675), each= 3)+0.15,
     rep(seq(0.3, 0.7, length.out=3)+0.06, 2),
     density = rep(c(40, NA), each= 3),
     col= c(rep(NA, 3), ex$col))
par(mar= c(2,2.6,0.25,0.1),
    cex.axis= 0.6, 
    cex.lab= 0.7,
    cex= 0.66)
for(mod in c("Add.", "Mult."))
{
  ex[, {
    y <- log2FoldChange
    z <- median_L
    if(mod=="Add.")
    {
      x <- 2^additive
      y <- 2^y
      z <- 2^z
      ylab <- "Activity"
      par(mgp= c(1.4, 0.5, 0))
    }else if(mod=="Mult.")
    {
      x <- multiplicative
      ylab <- "Multiplicative"
      par(mgp= c(1, 0.5, 0))
    }
    bar <- barplot(c(x, y), 
                   col= c(col, "lightgrey"),
                   xaxt= "n",
                   ylab= "Activity",
                   space= 1,
                   xlim= c(0.25,4))
    barplot(z, 
            col= "black",
            density= 40,
            add= T,
            space= 1)
    text(bar, 
         par("usr")[3]-strheight("M", cex= 0.5), 
         c(mod, "Obs."),
         pos= 2, 
         offset= -0.25, 
         xpd= T, 
         srt= 45, 
         cex= 0.8)
  }, .(L, R)]
}

# Boxplot residuals
par(mar= c(4,6,2,3),
    mgp= c(1.75, 0.5, 0),
    cex.axis= 1, 
    cex.lab= 1,
    lwd= 1)
vl_boxplot(dat[, .(log2FoldChange-additive,
                   log2FoldChange-multiplicative,
                   log2FoldChange-predicted)],
           ylab= "Residuals (log2)",
           xaxt= "n")
text(1:3, 
     par("usr")[3]-strheight("M", cex= 0.5), 
     c("Additive", "Multiplicative", "Linear model"),
     pos= 2, 
     offset= -0.25, 
     xpd= T, 
     srt= 45)
abline(h= 0, lty= 2)
dev.off()
