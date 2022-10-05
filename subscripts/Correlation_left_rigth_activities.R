setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
dat <- fread("Rdata/final_300bp_enhancer_features.txt")[, .(ID, group, dev_log2FC_TWIST)]
dat[, TWIST_class:= fcase(group=="control", "Control",
                          dev_log2FC_TWIST<=2 | is.na(dev_log2FC_TWIST), "Inactive",
                          dev_log2FC_TWIST<=6, "Weak",
                          dev_log2FC_TWIST<=8, "Medium",
                          dev_log2FC_TWIST>8, "Strong")]
dat[, TWIST_class:= factor(TWIST_class, 
                           c("Control", "Inactive", "Weak", "Medium", "Strong"))]
dat[, TWIST_col:= c("grey20", "#4D9221", "#B8E186", "#F1B6DA", "#C51B7D")[TWIST_class]]

dat <- dat[unique(lib[, .(ID= L, median_L)]), on= "ID", nomatch= NULL]
dat <- dat[unique(lib[, .(ID= R, median_R)]), on= "ID", nomatch= NULL]
dat <- dat[order(group=="control")]

#------------------------------------------------------#
# PLOT
#------------------------------------------------------#
pdf("pdf/draft/Correlation_left_rigth_activities.pdf", 
    width = 3, 
    height = 3)
# Scatter plot
par(mar= c(3.5, 3, 0.5, 0.5), 
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2,
    bty= "n")
dat[, {
  # Scatter plot
  plot(median_R,
       median_L,
       pch= 16,
       cex= 0.5,
       col= adjustcolor(TWIST_col, 0.7),
       las= 1,
       xlab= "3' individual activity (log2)",
       ylab= "5' individual activity (log2)")
  # Linear model
  .lm <- lm(median_L~median_R)
  abline(.lm, lty= "11")
  
  # Legend
  leg <- unique(.SD[, TWIST_col, keyby= TWIST_class])
  leg[, legend("topleft",
               legend = c(paste0("PCC= ", round(cor.test(median_L, median_R)$estimate, 2)),
                          as.character(rev(TWIST_class))),
               pch= c(NA, 
                      rep(16, length(TWIST_class))),
               col= c("black", 
                      rev(TWIST_col)),
               lty= c("11",
                      rep(NA, length(TWIST_class))),
               cex= 0.7,
               seg.len= 0.5,
               bty= "n",
               y.intersp= 0.8)]
}]
dev.off()



