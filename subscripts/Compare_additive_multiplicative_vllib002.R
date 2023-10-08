setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]

# Linear ----
model <- readRDS("db/linear_models/lm_vllib002.rds")
adj.rsqlm <- summary(model)$adj.r.square
eq <- model$coefficients
eq[c(1,4)] <- round(eq[c(1,4)], 2)
eq[2:3] <- round(eq[2:3], 1)
eq <- paste0("Predicted= ", eq[1], "+", eq[2], "*5'+", eq[3], "*3'", eq[4], "*5':3'")

# Melt data for plotting ----
setnames(dat,
         c("predicted", "log2FoldChange"),
           c("Linear model", "Observed"))
pl <- melt(dat, 
           id.vars = c("L", "R", "Observed", "actPair"), 
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl[, Rsq:= {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  1-(((1-rsq)*(.N-1))/(.N-2-1))
}, variable]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)",
                   "Linear model"= "Predicted activity (log2)"), variable]

# Check best predictor ----
t1 <- pl[(actPair) & variable!="Linear model", .(variable[which.min(abs(residuals))]), .(L, R)]$V1
t1 <- table(droplevels(t1))
t2 <- pl[(actPair), .(variable[which.min(abs(residuals))]), .(L, R)]$V1
t2 <- table(droplevels(t2))

# Clip extremes for ploting ----
clip <- c(0.005, 0.999)

# Plot ----
pdf("pdf/draft/Compare_add_mult_vllib002.pdf",
    width = 15,
    height = 3)
layout(matrix(1:5,
              nrow= 1))
par(mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
pl[, {
  # Model
  clip.x <- between(value,
                    quantile(value, clip[1]),
                    quantile(value, clip[2]), incbounds = T)
  clip.y <- between(Observed,
                    quantile(Observed, clip[1]),
                    quantile(Observed, clip[2]), incbounds = T) 
  .SD[clip.x & clip.y, {
    smoothScatter(value,
                  Observed,
                  xlab= xlab,
                  ylab= "Combined activity (log2)",
                  main= variable,
                  xaxt= "n",
                  yaxt= "n",
                  col= adjustcolor(blues9[9], .3))
  }]
  axis(1, padj = -1.25)
  axis(2)
  # Density active pairs
  z <- MASS::kde2d(value[(actPair)],
                   Observed[(actPair)], 
                   n = 50)
  contour(z,
          lwd= 0.25,
          drawlabels= FALSE,
          nlevels= 6,
          add= TRUE,
          xpd= NA)
  leg.x <- par("usr")[1]+strwidth("M")*0.75
  leg.y <- par("usr")[4]-strheight("M")*.75
  points(x = rep(leg.x, 2), # Density legend
         y = rep(leg.y, 2),
         lwd= 0.5,
         cex= c(0.6, 1.25))
  text(x = leg.x,
       y = leg.y,
       pos= 4,
       labels= "enh./enh. pairs",
       cex= 7/12)
  # R2
  vl_plot_coeff(value = Rsq,
                inset= c(-0.075, 0.06),
                cex= 7/12)
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
  lw <- diff(grconvertX(0:1, "line", "user"))
  lh <- diff(grconvertY(0:1, "line", "user"))
  clip(par("usr")[1]+lw,
       par("usr")[2]-lw,
       par("usr")[3]+lh,
       par("usr")[4]-lh)
  abline(0, 1, lty= "11")
  print(paste0(.GRP, "/", .NGRP))
}, .(variable, xlab, Rsq)]
# Performances
par(mai= c(.9,1.25,.9,1.25),
    pty= "m",
    mgp= c(1, .25, 0))
unique(pl[, .(variable, Rsq)])[, {
  bar <- barplot(Rsq,
                 ylab= "Adj. R2 (goodness of fit)",
                 col= "lightgrey")
  vl_tilt_xaxis(bar, 
                labels = variable)
  
}]
# Boxplot examples
par(mai= c(.9,1.25,.9,1.1),
    mgp= c(.75, .25, 0))
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
# Pie charts
plot.new()
par(mai= rep(1, 4),
    pty= "s")
pie(t1,
    labels = paste0(names(t1), " (", round(t1/sum(t1)*100), "%) n= ", formatC(t1, big.mark = ",")),
    cex= 8/12,
    lwd= .5,
    xpd= NA)
pie(t2,
    labels = paste0(names(t2), " (", round(t2/sum(t2)*100), "%), n= ", formatC(t2, big.mark = ",")),
    cex= 8/12,
    lwd= .5,
    xpd= NA)
dev.off()