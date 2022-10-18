setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
model <- readRDS("db/linear_models/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_lm.rds")
rsq <- round(summary(model)$r.squared, 2)

pdf("pdf/draft/CV_linear_model_vllib002.pdf",
    height = 3,
    width = 3)
par(las= 1,
    mar= c(3,3,1,0.5),
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
                ylab= "Predicted activity (log2)")
  abline(0,1,lty=2)
  text(par("usr")[1],
       par("usr")[4]-strheight("M"),
       paste0("R2= ", rsq),
       pos= 4)
}]
dev.off()