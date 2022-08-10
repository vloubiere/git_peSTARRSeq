setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

source("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/peAddFeatures_function.R")

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
dat <- fread("Rdata/final_300bp_enhancer_features.txt",
             sel= match(c("ID", "group", "dev_log2FC_TWIST"), peAddFeatures_list()))

dat <- dat[unique(lib[, .(ID= L, median_L)]), on= "ID", nomatch= NULL]
dat <- dat[unique(lib[, .(ID= R, median_R)]), on= "ID", nomatch= NULL]
dat[, TWIST_class:= fcase(group=="control", "Neg. ctl.",
                          dev_log2FC_TWIST<=2 | is.na(dev_log2FC_TWIST), "inact. (<=2/NA)",
                          dev_log2FC_TWIST<=6, "weak (2-6)",
                          dev_log2FC_TWIST<=8, "medium (6-8)",
                          dev_log2FC_TWIST>8, "strong (>8)")]
dat[, TWIST_class:= factor(TWIST_class, 
                           c("Neg. ctl.", "inact. (<=2/NA)", "weak (2-6)", "medium (6-8)", "strong (>8)"))]
dat[, TWIST_col:= c("grey20", "#4D9221", "#B8E186", "#F1B6DA", "#C51B7D")[TWIST_class]]
dat <- dat[order(group=="control")]

#------------------------------------------------------#
# PLOT
#------------------------------------------------------#
pdf("pdf/draft/Sup_Fig_1A.pdf", 
    width = 3, 
    height = 3)
# Scatter plot
par(mar= c(3.5, 3, 0.5, 0.5), 
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2)
dat[, {
  # Scatter plot
  plot(median_L, 
       median_R, 
       pch= 16,
       cex= 0.5,
       col= adjustcolor(TWIST_col, 0.7),
       las= 1,
       xlab= "5' individual activity (log2)",
       ylab= "3' individual activity (log2)",
       xaxt= "n",
       yaxt= "n")
  axis(1, 
       lwd= 0,
       lwd.ticks= 1)
  axis(2, 
       lwd= 0,
       lwd.ticks= 1)
  # Linear model
  abline(h= 0, lty= "11")
  abline(v= 0, lty= "11")
  .lm <- lm(median_R~median_L)
  abline(.lm)
  # Legend
  leg <- unique(.SD[, TWIST_col, keyby= TWIST_class])
  leg[, legend(par("usr")[1]+strwidth("M", cex= 0.1),
               par("usr")[4]-strheight("M", cex= 0.1),
               legend = TWIST_class,
               pch= 16,
               col= TWIST_col,
               box.lty= "11",
               cex= 0.5,
               seg.len= 1,
               bg= "white")]
  # PCC
  legend("bottomright",
         legend = paste0("R²= ", round(summary(.lm)$r.squared, 2),
                         " (PCC= ", round(cor.test(median_L, median_R)$estimate, 2), ")"),
         lty= 1,
         bty= "n",
         cex= 0.6,
         seg.len= 1)
}]
dev.off()



