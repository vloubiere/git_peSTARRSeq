setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
dat[, basal:= mean(log2FoldChange[ctlL & ctlR])]
dat[, additive:= log2(2^indL+2^indR-2^basal)]
dat[, multiplicative:= log2(2^indL*2^indR/2^basal)]

# Different models ----
## Additive
rsqAdd <- vl_model_eval(dat[!is.na(additive),log2FoldChange], dat[!is.na(additive),additive])$Rsquare
adj.rsqAdd <- 1-(((1-rsqAdd)*(nrow(dat)-1))/(nrow(dat)-2-1))

## Multiplicative
rsqMult <- vl_model_eval(dat$log2FoldChange, dat$multiplicative)$Rsquare
adj.rsqMult <- 1-(((1-rsqMult)*(nrow(dat)-1))/(nrow(dat)-2-1))

## Linear
# model <- readRDS("db/linear_models/lm_vllib002.rds")
model <- lm(log2FoldChange~indL*indR, dat)
dat[, predicted:= predict(model)]
dat[, residuals:= log2FoldChange-predicted]
adj.rsqlm <- summary(model)$adj.r.square
eq <- model$coefficients
eq[c(1,4)] <- round(eq[c(1,4)], 2)
eq[2:3] <- round(eq[2:3], 1)
eq <- paste0("Predicted= ", eq[1], "+", eq[2], "*5'+", eq[3], "*3'", eq[4], "*5':3'")

# Melt data for plotting ---
pl <- melt(dat, 
           id.vars = "log2FoldChange", 
           measure.vars = c("additive", "multiplicative", "predicted"))
pl[, residuals:= log2FoldChange-value]
pl[, Rsq:= switch(as.character(variable), 
                  "additive"= adj.rsqAdd,
                  "multiplicative"= adj.rsqMult,
                  "predicted"= adj.rsqlm), variable]
pl[, xlab:= switch(as.character(variable), 
                   "additive"= "5' + 3' activities (log2)",
                   "multiplicative"= "5' x 3' activities (log2)",
                   "predicted"= "Predicted activity (log2)"), variable]
pl[, title:= switch(as.character(variable), 
                       "additive"= "Additive model",
                       "multiplicative"= "Multiplicative model",
                       "predicted"= "Linear model"), variable]

pdf1 <- tempfile(fileext = ".pdf")
pdf2 <- tempfile(fileext = ".pdf")
  
# Barplot ----
pdf(pdf1,
    height = 1.875,
    width = 1)
par(las= 1,
    mar= c(4.5,3.25,1,0.75),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    cex= 8/12,
    cex.axis= 7/8,
    lwd= 0.75)
bar <- barplot(c(adj.rsqAdd, adj.rsqMult, adj.rsqlm),
               ylab= "Adj. R2 (goodness of fit)")
vl_tilt_xaxis(bar, 
              labels = c("Additive", "Multiplicative", "Linear model"))
dev.off()

# Scatterplots ----
lim <- c(-4, 14)
pdf(pdf2,
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
pl[, {
  # Model
  smoothScatter(value,
                log2FoldChange,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= xlab,
                ylab= "Combined activity (log2)",
                main= title,
                xlim= lim,
                ylim= lim,
                xaxt= "n",
                yaxt= "n",
                nrpoints = 0)
  axis(1, c(-4,0,6,12), c(-4,0,6,12))
  axis(2, c(-4,0,6,12), c(-4,0,6,12))
  legend('topleft',
         legend= bquote(Adj.~R^2 == .(round(Rsq, 2))),
         bty= "n",
         cex= 7/8)
  clip(-3.75,11,-3.75,11)
  abline(0, 1, lty= "11")
  if(variable=="predicted")
    text(mean(par("usr")[c(1,2)]),
         par("usr")[4],
         eq,
         pos= 3,
         xpd= T,
         offset= -0.1,
         cex= 0.7)
  # Diagnostic
  smoothScatter(value,
                residuals,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Fitted (log2)",
                ylab= "Residuals (log2)",
                main= title,
                ylim= c(-10, 6),
                xaxt= "n",
                nrpoints = 0)
  ticks <- if(variable=="additive")
    c(0,4,8) else if(variable=="multiplicative")
      c(0,6,12) else if(variable=="predicted")
        c(0,5,10)
  axis(1, ticks, ticks)
  abline(h= 0, lwd= 0.5)
  .l <- loess(residuals~value,
              .SD[sample(.N, 50000)][order(value)])
  lines(.l$x, .l$fitted, lty= "11")
  print(paste0(.GRP, "/", .NGRP))
}, .(variable, title, xlab, Rsq)]
dev.off()

pdftools::pdf_combine(c(pdf1, pdf2),
                      "pdf/draft/Compare_add_mult_vllib002.pdf")