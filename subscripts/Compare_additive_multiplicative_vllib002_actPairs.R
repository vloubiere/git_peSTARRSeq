setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import and compute basalMean
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat[, basalMean:= mean(log2FoldChange[ctlL & ctlR])]
actPairs <- readRDS("db/linear_models/FC_vllib002_actPairs_lm_predictions.rds")
dat <- dat[actPairs[, .(L, R)], on= c("L", "R")]

# Compute additive and multiplicative
dat[, `Additive model`:= log2(2^indL+2^indR-2^basalMean)]
dat[, `Multiplicative model`:= log2(2^indL*2^indR/2^basalMean)]

# Linear ----
model <- readRDS("db/linear_models/lm_vllib002_actPairs.rds")
adj.rsqlm <- summary(model)$adj.r.square
eq <- model$coefficients
eq[c(1,4)] <- round(eq[c(1,4)], 2)
eq[2:3] <- round(eq[2:3], 1)
eq <- paste0("Predicted= ", eq[1], "+", eq[2], "*5'+", eq[3], "*3'", eq[4], "*5':3'")

# Melt data for plotting ---
setnames(dat, "predicted", "Linear model")
pl <- melt(dat, 
           id.vars = "log2FoldChange", 
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= log2FoldChange-value]
pl[, Rsq:= {
  rsq <- vl_model_eval(log2FoldChange, value)$Rsquare
  1-(((1-rsq)*(.N-1))/(.N-2-1))
}, variable]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)",
                   "Linear model"= "Predicted activity (log2)"), variable]
# Clip extremes for ploting
clip <- c(0.005, 0.995)
pl <- pl[, {
  .SD[between(value,
              quantile(value, clip[1]),
              quantile(value, clip[2]), incbounds = T)
      & between(log2FoldChange,
                quantile(value, clip[1]),
                quantile(value, clip[2]), incbounds = T)]
}, variable]

pdf("pdf/draft/Compare_add_mult_vllib002_actPairs.pdf",
    width = 6.25,
    height = 2.25)
layout(matrix(1:4,
              nrow= 1),
       widths= c(0.5,1,1,1))
par(font.main= 1,
    las= 1,
    lend= 2,
    cex= 8/12,
    cex.axis= 7/8,
    cex.main= 1,
    mar= c(5,3,2,.5),
    mgp= c(1.5, 0.325, 0),
    tcl= -0.2,
    bty= "n",
    lwd= .5)
unique(pl[, .(variable, Rsq)])[, {
  bar <- barplot(Rsq,
                 ylab= "Adj. R2 (goodness of fit)",
                 lwd= 1)
  vl_tilt_xaxis(bar, 
                labels = variable)
  
}]
par(lwd= 1)
pl[, {
  # Model
  smoothScatter(value,
                log2FoldChange,
                xlab= xlab,
                ylab= "Combined activity (log2)",
                main= variable,
                xaxt= "n",
                yaxt= "n",
                col= adjustcolor(blues9[9], .3))
  axis(1)
  axis(2)
  # Legends
  vl_plot_R2(rsquare = Rsq,
             inset= c(-0.1, 0))
  lw <- diff(grconvertX(0:1, "line", "user"))
  lh <- diff(grconvertY(0:1, "line", "user"))
  clip(par("usr")[1]+lw,
       par("usr")[2]-lw,
       par("usr")[3]+lh,
       par("usr")[4]-lh)
  abline(0, 1, lty= "11")
  if(variable=="Linear model")
    text(mean(par("usr")[c(1,2)]),
         par("usr")[4],
         eq,
         pos= 3,
         xpd= T,
         offset= -0.1,
         cex= 0.7)
  print(paste0(.GRP, "/", .NGRP))
}, .(variable, xlab, Rsq)]
dev.off()