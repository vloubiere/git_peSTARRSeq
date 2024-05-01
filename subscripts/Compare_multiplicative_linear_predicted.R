setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data ----
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")

# Plot diagnostic ----
pdf("pdf/draft/Compare_mult_linear_predictions_vllib002.pdf", 6, 3)
vl_par(mfrow= c(1,2),
       mgp= c(1, .35, 0),
       font.main= 1,
       pty= "s")
dat[, {
  ylim <- range(c(log2FoldChange-`Multiplicative model`,
                  log2FoldChange-`Linear model`))
  smoothScatter(`Multiplicative model`,
                log2FoldChange-`Multiplicative model`,
                xlab= "Predicted value (log2)",
                ylab= "Residuals (log2)",
                main= "Multiplicative model",
                col= adjustcolor(blues9[9], .3),
                xaxt= "n")
  axis(1, padj= -1.25)
  ss <- smooth.spline(`Multiplicative model`,
                      log2FoldChange-`Multiplicative model`)
  abline(h= 0, lty= "13", lwd= .5)
  lines(ss, col= "red")
  smoothScatter(`Linear model`,
                log2FoldChange-`Linear model`,
                xlab= "Predicted value (log2)",
                ylab= "Residuals (log2)",
                main= "Multiplicative model\nwith interaction term",
                col= adjustcolor(blues9[9], .3),
                xaxt= "n")
  axis(1, padj= -1.25)
  ss <- smooth.spline(`Linear model`,
                      log2FoldChange-`Linear model`)
  abline(h= 0, lty= "13", lwd= .5)
  lines(ss, col= "red")
}]
dev.off()