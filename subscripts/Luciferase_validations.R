setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- merge(luc,
             lib, 
             by= c("L", "R"),
             suffixes= c("_luc", "_STARR"))
dat[, actCol:= adjustcolor(c("grey0", 
                             "grey25", 
                             "grey50", 
                             "grey75", 
                             "royalblue2", 
                             "cornflowerblue", 
                             "purple", 
                             "magenta", 
                             "#74C27A")[actClass], 0.5)]
dat <- dat[actClass %in% c("ctl./ctl.", "ctl./enh.", "enh./ctl.", "enh./enh.")]
dat[, actClass:= droplevels(actClass)]

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
       col= actCol,
       pch= 16,
       cex= 0.8)
  segments(log2FoldChange_STARR,
           log2FoldChange_luc-sd,
           log2FoldChange_STARR,
           log2FoldChange_luc+sd,
           col= actCol)
  .lm <- lm(log2FoldChange_luc~log2FoldChange_STARR, dat)
  abline(.lm, lty= "11")
  leg <- unique(dat[order(actClass), .(actClass, actCol)])
  leg[, legend("topleft",
               legend = c(paste0("R2= ", round(summary(.lm)$r.squared, 2)),
                          as.character(rev(actClass))),
               bty= "n",
               col= c("black", 
                      rev(actCol)),
               pch= c(NA, 
                      rep(16, length(actClass))),
               lty= c("11", 
                      rep(NA, length(actClass))),
               seg.len= 0.5,
               cex= 0.8)]
}]
dev.off()



