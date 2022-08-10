setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- merge(luc,
             lib, 
             by= c("L", "R"),
             suffixes= c("_luc", "_STARR"))
dat[, class_act:= droplevels(class_act)]
dat[, col_act:= adjustcolor(col_act, 0.5)]

pdf("pdf/draft/Figure_1FG.pdf", 
    height = 3, 
    width = 4.8)
layout(matrix(1:2, ncol= 2), widths = c(1, 0.65))
par(las= 1,
    mar= c(3.5, 2.75, 0.5, 0.5),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2)
dat[, {
  plot(log2FoldChange_STARR, 
       log2FoldChange_luc,
       xlab= "pe-STARR-Seq activity (log2)",
       ylab= "Normalized luc. activity (log2)",
       ylim= c(-0.8, 6.6),
       col= col_act,
       pch= 16,
       xaxt= "n",
       yaxt= "n")
  segments(log2FoldChange_STARR,
           log2FoldChange_luc-sd,
           log2FoldChange_STARR,
           log2FoldChange_luc+sd,
           col= col_act)
  axis(1, lwd= 0, lwd.ticks= 1)
  axis(2, lwd= 0, lwd.ticks= 1)
  leg <- unique(dat[order(class_act), .(class_act, col_act)])
  leg[, legend("topleft",
               legend = class_act,
               bty= "n",
               col= col_act,
               pch= 19,
               cex= 0.7)]
  .lm <- lm(log2FoldChange_luc~log2FoldChange_STARR, dat)
  abline(.lm, lty=2)
  legend("bottomright",
         legend = paste0("R²= ", round(summary(.lm)$r.squared, 2), 
                         " (PCC= ", round(cor.test(log2FoldChange_STARR, log2FoldChange_luc)$estimate, 2), ")"),
         bty= "n",
         lty= 2,
         cex= 0.7)
  par(mar= c(3.5, 3.5, 0.5, 0.5))
  vl_boxplot(log2FoldChange_luc~class_act,
             .SD,
             violin= T,
             compute_pval= list(c(1,2), c(2,4), c(3,4)),
             violcol= unique(col_act),
             ylab= "Normalized luc. activity (log2)",
             tilt.names= T, 
             pval_offset= 0.06,
             ylim= c(-0.5, 8))
}]
dev.off()
file.show("pdf/draft/Figure_1FG.pdf")
