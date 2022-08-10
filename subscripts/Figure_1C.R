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
dat[, class:= fcase(group=="control", "Random control sequences",
                    default = "Candidates")]
dat[, class:= factor(class, c("Random control sequences", "Candidates"))]
dat[, class_col:= c("grey20", "grey80")[class]]
dat <- dat[order(group=="control")]

#------------------------------------------------------#
# PLOT
#------------------------------------------------------#
pdf("pdf/draft/Figure_1C.pdf", 
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
       col= adjustcolor(class_col, 0.7),
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
  leg <- unique(.SD[, .(class, class_col)])
  leg[, legend(par("usr")[1]+strwidth("M", cex= 0.1),
               par("usr")[4]-strheight("M", cex= 0.1),
               legend = class,
               pch= 16,
               col= class_col,
               box.lty= "11",
               cex= 0.5,
               seg.len= 1,
               bg= "white",
               border= NA)]
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



