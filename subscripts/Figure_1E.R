setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- lib[vllib=="vllib002"]
dat <- dat[class %in% c("ctl./ctl.", "ctl./enh.", "enh./ctl.", "enh./enh.")]
dat[, class:= droplevels(class)]

pdf("pdf/draft/Figure_1E.pdf", width = 3, height = 4.5)
par(las= 2, mar= c(6,4,1,1))
vl_boxplot(log2FoldChange~class,
           dat,
           violin= T, 
           compute_pval = list(c(1,2), c(1,3), c(2,4), c(3,4)),
           violcol = dat[!is.na(col), col[1], keyby= class]$V1,
           ylab= "Activity (log2)", 
           ylab.line = 2, 
           wilcox.alternative = "less", 
           pval_adj= 0.06)
abline(h= 0, lty= 2)
dev.off()


