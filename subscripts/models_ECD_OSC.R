setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
ECD <- readRDS("db/linear_models/FC_DSCP_ECD_full_dataset_lm_predictions.rds")
ECD <- ECD[!is.na(group)] # Select ECD enhancer pairs and controls
OSC <- readRDS("db/linear_models/FC_DSCP_OSC_full_dataset_lm_predictions.rds")
OSC <- OSC[!is.na(group)] # Select OSC enhancer pairs and controls
dat <- list("S2 cells + ecdysone"= ECD,
            "OSC cells"= OSC)

# Import models ----
models <- list("S2 cells + ecdysone"= readRDS("db/linear_models/lm_DSCP_ECD_full_dataset.rds"),
               "OSC cells"= readRDS("db/linear_models/lm_DSCP_OSC_full_dataset.rds"))

# Fitted curves ----
mult <- data.table(indL= seq(0, 10, .1))
mult[, indR:= indL]
mult[, add:= log2(2^indL+2^indR-1)]
mult[, `S2 cells + ecdysone`:= predict(models$`S2 cells + ecdysone`, .SD)]
mult[, `OSC cells`:= predict(models$`OSC cells`, .SD)]

# Plot ----
Cc <- c("tomato", "royalblue", "limegreen", "grey30")

pdf("pdf/draft/super_additivity_OSC_ECD.pdf", 9, 3)
mat <- matrix(c(1,1,1,1,2,4,3,5,6,6,6,6), nrow= 2)
layout(mat)
vl_par(mgp= c(1, .35, 0),
       font.main= 1,
       cex.main= 8/12)
for(cdition in c("S2 cells + ecdysone", "OSC cells"))
{
  # random order before plotting ----
  set.seed(1)
  pl <- dat[[cdition]]
  pl <- pl[sample(seq(nrow(pl)), nrow(pl))]
  pl[, {
    # Scatter plot
    par(pty= "s",
        mai= rep(.9 ,4))
    plot(`Additive model`,
         log2FoldChange,
         pch= 16,
         col= adjustcolor(Cc[group], .4),
         cex= .2,
         main= cdition,
         xaxt= "n",
         xlab= "Predicted additive (log2)",
         ylab= "Combined activity (log2)")
    axis(1, padj= -1.25)
    legend("topleft",
           legend= levels(group),
           col= Cc,
           pch= 16,
           bty= "n",
           cex= 0.5,
           xpd= T)
    # Add fitted lines
    clip(min(`Additive model`),
         max(`Additive model`),
         min(log2FoldChange),
         max(log2FoldChange))
    lines(mult$add,
          mult[[cdition]])
    abline(0, 1, lty= "11")
    
    # Add Legends
    legend("bottomright", 
           lty= c(3, 1),
           legend= c("Additive", "Fit. mult."),
           bty= "n",
           cex= .4,
           xpd= T)
    
    # Individual plots
    xl <- par("usr")[1:2]
    yl <- par("usr")[3:4]
    par(mai= c(.45,.45,.45,.45),
        cex.axis= 5/12,
        lwd= .75)
    pl[, {
      plot(`Additive model`,
           log2FoldChange,
           pch= 16,
           col= adjustcolor(Cc[group], .2),
           cex= .2,
           main= group,
           xaxt= "n",
           xlab= "Predicted additive (log2)",
           ylab= "Combined activity (log2)",
           xlim= xl,
           ylim= yl)
      axis(1, padj= -2.75)
      # Add fitted lines
      clip(min(`Additive model`),
           max(`Additive model`),
           min(log2FoldChange),
           max(log2FoldChange))
      lines(mult$add,
            mult[[cdition]])
      abline(0, 1, lty= "11")
    }, group]
    
    # Boxplot quantif
    par(mai= c(.9 , 1.2, 1.1, 1.3),
        pty= "m",
        cex.axis= 7/12,
        lwd= 1)
    vl_boxplot(log2FoldChange-`Additive model`~group,
               .SD,
               tilt.names= T,
               show.empty= F,
               compute.pval= list(c(1,4)),
               col= adjustcolor(Cc, .5),
               ylab= "Residuals (log2)\nobserved-pred. additive",
               lwd= .75)
    abline(h= 0, lty= "11")
    .SD
  }]
}
dev.off()