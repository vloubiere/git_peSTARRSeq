setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import and select active pairs
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")

# Melt data for plotting ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat, 
           id.vars = c("L", "R", "Observed", "actPair"), 
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl <- pl[, {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  rsqAct <- vl_model_eval(Observed[(actPair)], value[(actPair)])$Rsquare
  .(Rsq= 1-(((1-rsq)*(.N-1))/(.N-2-1)),
    RsqAct= 1-(((1-rsqAct)*(sum(actPair)-1))/(sum(actPair)-2-1)))
}, variable]


# Plot ----
pdf("pdf/draft/Compare_rsq_models.pdf",
    width = 6,
    height = 3)
par(mfrow= c(1,2),
    mai= c(.9,1.25,.9,1.25), 
    mgp= c(1.25, .25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2)
unique(pl[, .(variable, Rsq, RsqAct)])[, {
  bar <- barplot(Rsq,
                 ylab= "Adj. R2 (goodness of fit)",
                 col= "lightgrey")
  vl_tilt_xaxis(bar, 
                labels = variable)
  bar <- barplot(RsqAct,
                 ylab= "Adj. R2 enh./enh. pairs",
                 col= "#FF7F00")
  vl_tilt_xaxis(bar, 
                labels = variable)
}]
dev.off()