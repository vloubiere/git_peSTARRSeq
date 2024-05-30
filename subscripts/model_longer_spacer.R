setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/linear_models/FC_DSCP_long_spacer_act_pairs_lm_predictions.rds")
  
# Import models ----
model <- readRDS("db/linear_models/lm_DSCP_long_spacer_act_pairs.rds")

# Fitted curves ----
mult <- data.table(indL= seq(0, 10, .1))
mult[, indR:= indL]
mult[, add:= log2(2^indL+2^indR-1)]
mult[, `linear model`:= predict(model, .SD)]

# Plot ----
pdf("pdf/draft/model_long_spacer.pdf", 4, 3)
layout(matrix(1:2, ncol= 2), widths = c(3, 1))
vl_par(mgp= c(1, .35, 0),
       font.main= 1,
       pty= "s",
       mai= rep(.9 ,4))
dat[, {
  # Scatter plot
  smoothScatter(`Additive model`,
                log2FoldChange,
                xlab= "Predicted additive (log2)",
                ylab= "Combined activity (log2)",
                xaxt= "n",
                col= adjustcolor(blues9[9], .3))
  axis(1,
       padj = -1.25,
       gap.axis = 0)
  # Fitting lines
  clip(min(`Additive model`, na.rm= T),
       max(`Additive model`, na.rm= T),
       min(log2FoldChange),
       max(log2FoldChange))
  lines(mult$add,
        mult$`linear model`,
        col= "red")
  abline(0,
         1,
         lty= "11")
  # Add Legends
  legend("bottomright", 
         lty= c(3, 1),
         col= c("black", "red"),
         legend= c("Additive", "Fit. mult."),
         bty= "n",
         cex= .4,
         xpd= T)
  # Residuals ----
  par(mai= c(.9, .5, .9, .1),
      pty= "m")
  vl_boxplot(log2FoldChange-`Additive model`,
             tilt.names = T,
             ylab= "Residuals (log2)\nobserved-pred. additive",
             lwd= .75)
  abline(h= 0, lty= "13")
}]
dev.off()