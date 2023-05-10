setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")

# Select a weak left enhancer for which residuals (among all pairs with Right enhancers) are centered on 0
sel <- dat[!grepl("^control", L) & !grepl("^control", R) & between(indL, 2, 4.5) & padjL<0.05, 
           abs(diff(c(0, median(residuals)))), L][which.min(V1), L]
ex <- dat[L==sel]
# Select examples for which right enhancers all have similar activities and similar to left enhancer
pts <- ex[between(indR, indL[1]-0.5, indL[1]+0.5)]
# AND that are either higher or lower than expected
pts <- pts[c(which.min(abs(-1-residuals)),
             which.min(abs(0-residuals)),
             which.min(abs(1-residuals)))]
Cc <- c("cornflowerblue", "limegreen", "tomato")

pdf("pdf/draft/lm_residuals_examples_vllib002.pdf", 5, 3)
par(mar= c(3,5,1,2),
    mfrow= c(1,2),
    tcl= -0.2,
    las= 1,
    lwd= 0.5,
    mgp= c(2.5,0.5,0))
vl_boxplot(ex$residuals,
           ylab= "Observed/Expected (log2)", 
           violin = T,
           col= "lightgrey")
points(rep(1, 3), pts[, residuals], col= Cc, pch= 19)
par(tcl= -0.2,
    las= 1,
    lwd= 0.5,
    mgp= c(1.25,0.25,0))
bar <- barplot(c(pts$indL[1], 
          pts$indR[1],
          pts$predicted[1],
          pts$log2FoldChange[1],
          pts$indR[2],
          pts$predicted[2],
          pts$log2FoldChange[2],
          pts$indR[3],
          pts$predicted[3],
          pts$log2FoldChange[3]),
        space= rep(c(0.25,1.5,0.25), length.out= 10),
        col= c("lightgrey", rep(Cc, each= 3)),
        border= NA,
        names.arg= c("5' candidate",
                     rep(c("3' candidate",
                           "5'/3' expected",
                           "5'/3' observed"), 3)),
        horiz= T,
        xaxt= "n")
axis(1, c(0, 5), c(0, 5))
title(xlab= "Activity (log2)")
adj <- strwidth("M")*0.25
pts[, {
  lines(c(log2FoldChange+adj, 
          rep(max(c(log2FoldChange,predicted))+adj*2, 2),
          predicted+adj),
        rep(bar[c(c(4, 7, 10)[.GRP], c(3, 6, 9)[.GRP]), 1], each= 2),
        xpd= T)
}, .(log2FoldChange, predicted)]
dev.off()

file.show("pdf/draft/lm_residuals_examples_vllib002.pdf")
