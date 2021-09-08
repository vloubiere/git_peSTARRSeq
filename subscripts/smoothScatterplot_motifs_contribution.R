setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

if(!exists("dat"))
{
  dat <- data.table(file= c("db/linear_models/global_models/vllib002_dev_x_dev__predict_log2FoldChange__actL*R.rds",
                            "db/linear_models/global_models/vllib002_dev_x_dev__predict_log2FoldChange__actL*R+motL*R.rds",
                            "db/linear_models/global_models/vllib002_dev_x_dev__predict_log2FoldChange__motL*R.rds",
                            "db/linear_models/global_models/vllib016_hk_x_hk__predict_log2FoldChange__actL*R.rds",
                            "db/linear_models/global_models/vllib016_hk_x_hk__predict_log2FoldChange__actL*R+motL*R.rds",
                            "db/linear_models/global_models/vllib016_hk_x_hk__predict_log2FoldChange__motL*R.rds"))
  dat[, form:= rep(c("log2Foldchange= Left act. * Right act.", 
                     "log2Foldchange= Left act. * Right act. + Left mot. * Right mot. [...]",
                     "log2Foldchange= Left mot. * Right mot. [...]"), 2)]
  dat[, title:= rep(c("Developmental x Developmental", 
                      "Housekeeping x Housekeeping"), each= 3)]
  dat <- dat[, {
    .c <- readRDS(file)
    data <- .c$model[, 1] 
    predict <- predict(.c)
    pcc <- cor.test(data, predict)
    pval <- paste0(round(pcc$estimate, 2), 
                   cut(pcc$p.value, c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
                       c("****", "***", "**", "*", "N.S")))
    .(data, 
      predict, 
      pval)
  }, (dat)]
}

pdf("pdf/modeling/smoothScatter_motifs_contribution.pdf")
par(mfrow= c(2,3),
    las= 1,
    pty= "s")
dat[, {
  smoothScatter(data,
                predict,
                xlab= "log2FoldChange",
                ylab= "Prediction")
  title(title, 
        font.main= 1)
  title(form, 
        font.main= 1, 
        cex.main= 0.7, 
        line = 0.75)
  abline(0, 1)
  legend("topleft", 
         legend = paste0("PCC= ", pval), 
         bty= "n", 
         cex= 0.8)
  print("")
}, .(file, form, title, pval)]
dev.off()
