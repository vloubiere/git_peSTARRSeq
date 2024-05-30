setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import and select active pairs
dat <- readRDS("db/linear_models/FC_RpS12_focused_lm_predictions.rds")
dat <- dat[grepl("^hk", L) & grepl("^hk", R)]

# Melt data for plotting ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat, 
           id.vars = c("L", "R", "Observed", "actPair"), 
           measure.vars = c("Additive model", "Multiplicative model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)"), variable]

# Compute ajusted Rsquared ----
pl[, nPred:= switch(as.character(variable), 
                    "Additive model"= 2,
                    "Multiplicative model"= 2), variable]
pl <- pl[, Rsq:= {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  .(1-(((1-rsq)*(.N-1))/(.N-nPred-1)))
}, .(variable, nPred)]

# Plot ----
pdf("pdf/draft/Modelling_obs_vs_expected_focused_RpS12_lib.pdf",
    width = 6,
    height = 3)
par(mfrow= c(1,2),
    mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
# Scatterplots ----
pl[, {
  xlim <- range(value)
  if(xlim[1]<(-5))
    xlim[1] <- -5
  smoothScatter(value,
                Observed,
                xlab= xlab,
                ylab= "Combined activity (log2)",
                main= variable,
                xaxt= "n",
                yaxt= "n",
                xlim= xlim,
                col= adjustcolor(blues9[9], .3))
  axis(1,
       padj = -1.25,
       gap.axis = 0)
  axis(2)
  # Density active pairs
  z <- MASS::kde2d(value[(actPair)],
                   Observed[(actPair)],
                   n = 50)
  contour(z,
          lwd= 0.5,
          drawlabels= FALSE,
          col= "tomato",
          nlevels= 10,
          add= TRUE,
          xpd= NA)
  leg.x <- par("usr")[1]+strwidth("M")*0.75
  leg.y <- par("usr")[4]-strheight("M")*.75
  leg.s <- c(.1, .33, .66, 1)
  points(x = rep(leg.x, length(leg.s)), # Density legend
         y = rep(leg.y, length(leg.s)),
         lwd= 0.5,
         cex= leg.s,
         col= "tomato")
  text(x = leg.x,
       y = leg.y,
       pos= 4,
       labels= "Enh./Enh. pairs",
       cex= 6/12,
       offset= .3,
       col= "tomato")
  # R2
  vl_plot_coeff(value = Rsq,
                inset= c(-0.075, 0.07),
                cex= 6/12)
  # Abline
  rg <- seq(-1, ifelse(variable=="Additive model", 7, 9))
  lines(rg, rg, lty= "11", col= "tomato")
  print(paste0(.GRP, "/", .NGRP))
}, .(variable, xlab, Rsq)]
dev.off()