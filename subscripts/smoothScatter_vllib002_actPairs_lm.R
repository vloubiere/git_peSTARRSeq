setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

dat <- readRDS("db/linear_models/FC_vllib002_actPAirs_lm_predictions.rds")
  
model <- readRDS("db/linear_models/lm_vllib002_actPAirs.rds")
Rsq <- summary(model)$adj.r.squared
eq <- model$coefficients
eq[c(1,4)] <- round(eq[c(1,4)], 2)
eq[2:3] <- round(eq[2:3], 1)
eq <- paste0("Predicted= ", eq[1], "+", eq[2], "*5'+", eq[3], "*3'", eq[4], "*5':3'")
title <- "Non-saturating active pairs"

#-----------------------------------------------#
# Scatterplots
#-----------------------------------------------#
pdf("pdf/draft/SmoothScatter_vllib002_actPairs_lm.pdf",
    height = 3.75,
    width = 6)
par(font.main= 1,
    las= 1,
    lend= 2,
    cex= 8/12,
    cex.axis= 7/8,
    cex.main= 1,
    mfcol= c(2,3),
    mar= c(4.1,4.1,2,2),
    mgp= c(1.5, 0.325, 0),
    tcl= -0.2,
    bty= "n")
dat[, {
  # Model
  smoothScatter(log2FoldChange, 
                predicted,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Combined activity (log2)",
                ylab= "Predicted activity (log2)",
                main= title,
                xaxt= "n",
                yaxt= "n")
  axis(1, c(0,4,8), c(0,4,8))
  axis(2, c(3,5,7), c(3,5,7))
  abline(0, 1, lty= "11")
  legend('topleft', 
         legend= bquote(Adj.~R^2 == .(round(Rsq, 2))),
         bty= "n",
         cex= 7/8)
  text(mean(par("usr")[c(1,2)]),
       par("usr")[4],
       eq,
       pos= 3,
       xpd= T,
       offset= -0.1,
       cex= 0.7)
  # Diagnostic
  smoothScatter(log2FoldChange,
                residuals,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Fitted (log2)",
                ylab= "Residuals (log2)",
                main= title)
  abline(h= 0, lty= "11")
  .l <- loess(residuals~predicted,
              .SD[sample(.N, 50000)][order(predicted)])
  lines(.l$x, .l$fitted, col= "red")
  print(paste0(.GRP, "/", .NGRP))
}]
dev.off()

file.show("pdf/draft/SmoothScatter_vllib002_actPairs_lm.pdf")
