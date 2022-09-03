setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
options(scipen= 1)
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
dat <- model$pred[dat, on= c("L", "R"), nomatch= NULL]
# Select examples
ex <- dat[list("dev_strong_B_00302",
               c("dev_medium_C_00480", "dev_medium_C_00466", "dev_medium_B_00585")), on= c("L", "R")]
ex[, col:= c("tomato", "limegreen", "cornflowerblue")]

#-----------------------------------------------#
# Plot
#-----------------------------------------------#
Cc <- colorRampPalette(c("white", "grey90", "grey60", "grey20"))
pdf("pdf/draft/Figure_2A.pdf",
    height = 3,
    width = 12)
layout(matrix(c(1,2,3,4,9,
                1,2,3,5,9,
                1,2,3,6,9,
                1,2,3,7,9,
                1,2,3,8,9), 
              byrow = T,
              ncol= 5),
       widths = c(1,1,1,0.325,0.5),
       heights =  c(0.06,0.2,0.2,0.2,0.02))
par(las= 1,
    mar= c(2.35,3,1.4,0.5),
    mgp= c(1.25, 0.25, 0),
    tcl= -0.2,
    cex= 1)
# Additive
smoothScatter(dat[, .(`Expected additive (log2)`= additive, 
                      `Activity (log2)`= log2FoldChange)],
              main= "Additive",
              colramp = Cc)
abline(0,1,lty=2)
legend(par("usr")[1]-strwidth("M", cex= 0.6),
       par("usr")[4],
       c(model$eq, 
         paste0("R²= ", round(vl_model_eval(dat$log2FoldChange, 
                                            dat$additive)$Rsquare, 2))),
       bty= "n",
       cex= 1)
# Multiplicative
smoothScatter(dat[, .(`Expected multiplicative (log2)`= multiplicative, 
                      `Activity (log2)`=log2FoldChange)],
              main= "Multiplicative",
              colramp = Cc)
abline(0,1,lty=2)
legend(par("usr")[1]-strwidth("M", cex= 0.6),
       par("usr")[4],
       c(model$eq, 
         paste0("R²= ", round(vl_model_eval(dat$log2FoldChange, 
                                            dat$multiplicative)$Rsquare, 2))),
       bty= "n",
       cex= 1)
# lm
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
         paste0("R²= ", round(model$model$R2, 2))),
       bty= "n",
       cex= 1)

# Examples barplot
par(mar= c(0,0,0,0),
    lwd= 0.1)
plot.new()
segments(0.35, 
         seq(0.3, 0.7, length.out=3), 
         0.95,
         seq(0.3, 0.7, length.out=3))
rect(rep(c(0.425, 0.725), each= 3),
     rep(seq(0.3, 0.7, length.out= 3)-0.06, 2),
     rep(c(0.425, 0.725), each= 3)+0.15,
     rep(seq(0.3, 0.7, length.out=3)+0.06, 2),
     col= c(rep("gold", 3), ex$col))
par(mar= c(3,2,0.1,0.1),
    cex.axis= 0.6,
    mgp= c(1, 0.5, 0),
    cex.lab= 0.7,
    cex= 0.66)
ex[, {
  bar <- barplot(c(log2FoldChange,
                   median_L, 
                   median_R,
                   additive,
                   multiplicative,
                   predicted), 
                 col= c("black", "gold", col, grey.colors(3, start= 0.9, end= 0.4)),
                 xaxt= "n",
                 ylab= "Activity",
                 space= 1)
  segments(x0 = mean(bar[1:2]),
           y0 = log2FoldChange,
           x1 = par("usr")[2],
           y1 = log2FoldChange,
           lty= "99")
  text(bar,
       par("usr")[3]-strheight("M", cex= 0.5),
       c("Obs.", "5'", "3'", "Exp. Add.", "Exp. Mult.", "Exp. lm"),
       pos= 2,
       offset= -0.25,
       xpd= T,
       srt= 45,
       cex= 0.8)
}, .(L, R)]

# Boxplot residuals
par(mar= c(0,0,0,0),
    lwd= 0.1)
plot.new()
par(mar= c(6.5,5,1,0.5),
    mgp= c(2.25, 0.5, 0),
    cex.axis= 1/0.66, 
    cex.lab= 1/0.66,
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
     srt= 45,
     cex= 1/0.66)
abline(h= 0, lty= 2)
dev.off()
