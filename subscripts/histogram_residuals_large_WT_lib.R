setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")

# plot ----
pdf("pdf/draft/histogram_residuals.pdf", 
    width = 3.1,
    height = 2.75)
par(mai= c(0.75,0.75,0.75,0.75), 
    mgp= c(1, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2)
dat[,{
  plot.new()
  plot.window(xlim= c(-4.5,4.5),
              ylim= c(0,.5))
  text(0,
       .5,
       pos= 1,
       labels = "Super-additive",
       cex= 7/12,
       offset= 0.1)
  rect(par("usr")[1],
       0,
       -2,
       .5,
       col= adjustcolor("cornflowerblue", .3),
       border= NA)
  text(mean(c(par("usr")[1], -2)),
       .5,
       pos= 1,
       labels = "Weaker",
       cex= 7/12,
       offset= 0.1)
  rect(2,
       0,
       par("usr")[2],
       .5,
       col= adjustcolor("tomato", .3),
       border= NA)
  text(mean(c(2, par("usr")[2])),
       .5,
       pos= 1,
       labels = "Stronger",
       cex= 7/12,
       offset= 0.1)
  par(lwd= 0.5)
  hist(residuals,
       xlim= c(-4.5,4.5),
       freq = F,
       add= T)
  axis(1, padj = -1.25)
  axis(2)
  title(xlab= "Residuals (log2)\n(observed-predicted)")
  title(ylab= "Density")
}]
dev.off()