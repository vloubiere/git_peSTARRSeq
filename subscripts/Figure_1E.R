setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- lib[vllib=="vllib002"]
dat <- dat[class %in% c("ctl./ctl.", "enh./ctl.", "ctl./enh.", "enh./enh.")]
dat[, class:= droplevels(class)]

pdf("pdf/draft/Figure_1E.pdf", 
    width = 2, 
    height = 3)
par(las= 2, 
    mar= c(3.5,3,0.5,2),
    tcl= -0.2,
    mgp= c(1.5, 0.5, 0))
vl_boxplot(log2FoldChange~class,
           dat,
           violin= T, 
           compute_pval = list(c(1,2), c(1,3), c(2,4), c(3,4)),
           violcol = dat[!is.na(col), col[1], keyby= class]$V1,
           ylab= "Activity (log2)", 
           ylab.line = 2, 
           wilcox.alternative = "less", 
           pval_offset= 0.05, 
           tilt.names= T,
           ylim= c(-5, 15))
abline(h= 0, 
       lty= 2)
dev.off()