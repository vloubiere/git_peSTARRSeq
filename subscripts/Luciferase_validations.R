setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- merge(luc,
             lib, 
             by= c("L", "R"),
             suffixes= c("_luc", "_STARR"))
dat[, class_act:= droplevels(class_act)]
dat[, col_act:= adjustcolor(col_act, 0.5)]

pdf("pdf/draft/Luciferase_validations.pdf", 
    height = 3, 
    width = 3)
par(las= 1,
    mar= c(3.5, 2.75, 0.5, 0.5),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    bty= "n")
dat[, {
  plot(log2FoldChange_STARR, 
       log2FoldChange_luc,
       xlab= "pe-STARR-Seq activity (log2)",
       ylab= "Norm. luciferase activity (log2)",
       ylim= c(-0.8, 6.6),
       col= col_act,
       pch= 16,
       cex= 0.8)
  segments(log2FoldChange_STARR,
           log2FoldChange_luc-sd,
           log2FoldChange_STARR,
           log2FoldChange_luc+sd,
           col= col_act)
  .lm <- lm(log2FoldChange_luc~log2FoldChange_STARR, dat)
  abline(.lm, lty= "11")
  leg <- unique(dat[order(class_act), .(class_act, col_act)])
  leg[, legend("topleft",
               legend = c(paste0("PCC= ", round(cor.test(log2FoldChange_STARR, log2FoldChange_luc)$estimate, 2)),
                          as.character(rev(class_act))),
               bty= "n",
               col= c("black", 
                      rev(col_act)),
               pch= c(NA, 
                      rep(16, length(class_act))),
               lty= c("11", 
                      rep(NA, length(class_act))),
               seg.len= 0.5,
               cex= 0.8)]
}]
dev.off()



