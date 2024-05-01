setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import and select active pairs
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")

# Linear ----
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
adj.rsqlm <- summary(model)$adj.r.square
eq <- model$coefficients
eq <- round(eq, 1)
eq <- paste0("Predicted= ", eq[1], "+", eq[2], "*5'+", eq[3], "*3'", eq[4], "*5':3'")

# Melt data for plotting ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat, 
           id.vars = c("L", "R", "Observed", "actPair"), 
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)",
                   "Linear model"= "Predicted activity (log2)"), variable]

# Compute ajusted Rsquared ----
pl[, nPred:= switch(as.character(variable), 
                   "Additive model"= 2,
                   "Multiplicative model"= 2,
                   "Linear model"= 3), variable]
pl <- pl[, Rsq:= {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  .(1-(((1-rsq)*(.N-1))/(.N-nPred-1)))
}, .(variable, nPred)]

# Plot ----
pdf("pdf/draft/Modelling_obs_vs_expected_large_WT_lib.pdf",
    width = 9,
    height = 3)
par(mfrow= c(1,3),
    mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
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
  # Equation
  if(variable=="Linear model")
    text(mean(par("usr")[c(1,2)]),
         par("usr")[4],
         eq,
         pos= 3,
         xpd= T,
         offset= .5,
         cex= 0.7)
  # Abline
  rg <- seq(-1, ifelse(variable=="Additive model", 7, 9))
  lines(rg, rg, lty= "11", col= "tomato")
  print(paste0(.GRP, "/", .NGRP))
}, .(variable, xlab, Rsq)]

# Boxplots residuals ----
par(mai= c(.9,1.25,.9,1.25),
    pty= "m",
    mgp= c(1, .25, 0))
pl[, {
  vl_boxplot(residuals,
             residuals[(actPair)],
             col= c(blues9[4], "tomato"),
             tilt.names= T,
             ylab= "Residuals\n(Observed-expected)",
             main= variable,
             violin= T,
             ylim= c(-2.5, 3.5),
             names= c("All", "Enh./Enh. pair"))
  abline(h= 0,
         lty="11")
  .SD
}, variable]

# Boxplot examples ----
dat[(actPair), {
  box <- vl_boxplot(.SD,
                    notch = T,
                    # compute_pval= list(c(1,2), c(3,4)),
                    tilt.names= T,
                    ylab= "Activity (log2)")$stats
  lim <- par("usr")[c(3,4)]
  sub <- .SD[apply(.SD, 1, function(x) all(x>lim[1] & x<lim[2]))
             & between(Observed,
                       box[1,4],
                       box[5,4],
                       incbounds = T)]
  sel <- seq(min(sub$Observed),
             max(sub$Observed),
             length.out= 25)
  sub <- lapply(sel, function(x) sub[which.min(abs(Observed-x))])
  sub <- rbindlist(sub)
  sub[, {
    lines(1:4,
          unlist(.BY),
          col= adjustcolor("grey10", .5),
          lwd= 0.5)
    points(1:4,
           unlist(.BY),
           col= adjustcolor("grey10", .5),
           pch= 16,
           cex= 0.4)
  }, (sub)]
  .SD
}, .SDcols= c("Additive model", "Multiplicative model", "Linear model", "Observed")]
dev.off()