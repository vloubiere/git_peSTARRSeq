setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

models <- readRDS("db/linear_models/lasso_vllib002_residuals.rds")
dat <- readRDS("db/linear_models/FC_lasso_vllib002_residuals_predictions.rds")

pdf("pdf/draft/residuals_prediction_lasso.pdf",
    height = 3, 
    width = 6)
par(mfrow= c(1,2),
    mai= rep(.9, 4), 
    mgp= c(1, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
for(var in c("log2FoldChange", "residuals"))
{
  # barplot(models[[var]]$CV_rsqs$rsq, 
  #         ylab= "R2 (goodness of fit)", 
  #         ylim= c(0, .25),
  #         xlab= "Cross-validations")
  # title(main= var)
  x <- dat[[var]]
  y <- dat[[switch(var,
                   "log2FoldChange"= "predicted_lasso",
                   "residuals"= "predicted_residuals_lasso")]]
  smoothScatter(x, 
                y,
                xlab= "Observed (log2)",
                ylab= "Predicted (LASSO)",
                col= adjustcolor(blues9[9], .3),
                main= var,
                xaxt= "n")
  axis(1, padj= -1.25)
  vl_plot_coeff(type = "rsq",
                value = vl_model_eval(x, y)$Rsquare,
                cex= 7/12,
                inset= c(-0.1, 0))
}
dev.off()