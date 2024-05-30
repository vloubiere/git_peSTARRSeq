setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("Rdata/homotypic_validations_luciferase_final_table.rds")

# Plot ----
Cc <- "grey50"
lims <- c(3, 7.5)
pdf("pdf/draft/luciferase_homotypic_pairs.pdf", 3, 3)
vl_par(mgp= c(.75, .25, 0))
dat[, {
  # Scatter plot
  par(pty= "s",
      mai= rep(.9 ,4))
  plot(`Additive model`,
       log2FoldChange,
       pch= 16,
       col= Cc,
       xaxt= "n",
       xlab= "Predicted additive (log2)",
       ylab= "Combined activity (log2)",
       cex= .7,
       xlim= lims,
       ylim= lims)
  arrows(`Additive model`,
         log2FoldChange,
         `Additive model`,
         log2FoldChange+sd,
         angle = 90,
         length = .0125,
         lwd= .5)
  arrows(`Additive model`,
         log2FoldChange,
         `Additive model`,
         log2FoldChange-sd,
         angle = 90,
         length = .0125,
         lwd= .5)
  axis(1, padj= -1.25)
  clip(3, 7, 3, 7)
  abline(0, 1)
  .SD
}]
dev.off()