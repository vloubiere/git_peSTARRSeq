setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")
model <- readRDS("db/linear_models/lm_vllib002.rds")
rsq <- round(summary(model)$r.squared, 2)
eq <- vl_model_equation(model, digits= 1)
eq <- gsub("log2FoldChange", "Act.", eq)
eq <- gsub("indL", "5'", eq)
eq <- gsub("indR", "3'", eq)
eq <- gsub("  ", " ", eq)
eq <- gsub(" * ", "*", eq, fixed = T)
eq <- gsub(" =", "=", eq, fixed = T)

pdf("pdf/draft/CV_linear_model_vllib002.pdf",
    height = 3.1,
    width = 3)
par(las= 1,
    mar= c(3,3,1.5,0.5),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    cex= 1,
    bty= "n")
# Additive
dat[, {
  smoothScatter(log2FoldChange,
                predicted,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Observed activity (log2)",
                ylab= "Predicted activity (log2)",
                xaxt= "n",
                yaxt= "n",
                xlim= c(-3.5, 10),
                ylim= c(-3.5, 10))
  # lines(loess(predicted~log2FoldChange))
  axis(1, at= c(-3,0,3,6,9), labels = c(-3,0,3,6,9))
  axis(2, at= c(-3,0,3,6,9), labels = c(-3,0,3,6,9))
  mtext(eq, line= 0.5, adj = 1, cex= 0.8)
  abline(0,1,lty=2)
  text(par("usr")[1],
       par("usr")[4]-strheight("M"),
       paste0("R2= ", rsq),
       pos= 4)
}]
dev.off()

file.show("pdf/draft/CV_linear_model_vllib002.pdf")
