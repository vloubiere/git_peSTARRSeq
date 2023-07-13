setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

models <- readRDS("db/linear_models/lasso_vllib002_residuals.rds")
dat <- readRDS("db/linear_models/FC_lasso_vllib002_residuals_predictions.rds")

pdf("pdf/draft/residuals_prediction_lasso.pdf",
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
for(var in c("log2FoldChange", "residuals"))
{
  barplot(models[[var]]$CV_rsqs$rsq, 
          ylab= "R2 (goodness of fit)", 
          ylim= c(0, .25),
          xlab= "Cross-validations")
  title(main= var)
  x <- dat[[var]]
  y <- dat[[switch(var,
                   "log2FoldChange"= "predicted_lasso",
                   "residuals"= "predicted_residuals_lasso")]]
  smoothScatter(x, 
                y,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Observed (log2)",
                ylab= "Predicted (lasso)")
  legend('topleft',
         legend = paste0("R2= ", round(vl_model_eval(x, y), 2)$Rsquare),
         bty= "n",
         cex= 7/8)
}
dev.off()