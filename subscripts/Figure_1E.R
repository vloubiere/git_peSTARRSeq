setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
dat <- lib[class_act %in% c("ctl./ctl.", "enh./ctl.", "ctl./enh.", "enh./enh.")]
dat[, class_act:= droplevels(class_act)]
setorderv(dat, "class_act")

pdf("pdf/draft/Figure_1E.pdf", 
    width = 2, 
    height = 3)
par(las= 2, 
    mar= c(3.5,3,0.5,2),
    tcl= -0.2,
    mgp= c(1.5, 0.5, 0),
    lty= 1)
vl_boxplot(log2FoldChange~class_act,
           dat,
           compute_pval = list(c(1,2), c(1,3), c(2,4), c(3,4)),
           col = adjustcolor(unique(dat$col_act), 0.5),
           ylab= "Activity (log2)",
           tilt.names= T,
           notch= T)
abline(h= 0, 
       lty= 2)
dev.off()