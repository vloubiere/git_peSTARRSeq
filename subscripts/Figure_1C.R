setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- lib[vllib=="vllib002"]
dat <- feat$add_feature(dat, feat$lib)
dat[, class_L := factor(ifelse(group_L=="control", "control", class_L), c("control", "inactive", "active"))]
dat[, class_R := factor(ifelse(group_R=="control", "control", class_R), c("control", "inactive", "active"))]

par(las= 2)
layout(matrix(c(c(1,4,3,2)), ncol=2), widths = c(0.25,1), heights = c(1,0.25))
vl_boxplot(median_L~class, 
           unique(dat[, .(L, 
                          median_L, 
                          class_L)]),
           horizontal= T, 
           violin= T)
abline(v= 0, lty= 2)
par(las= 1)
vl_boxplot(median_R~class, 
           unique(dat[, .(R, 
                          median_R, 
                          class= factor(ifelse(group_R=="control", "control", class_R),
                                        c("control", "inactive", "active")))]),
           violin= T)
abline(h= 0, lty= 2)



pdf("pdf/draft/Figure_1C.pdf", width = 2.5, height = 4.5)
par(las= 2, mar= c(6,4,1,1))
vl_boxplot(list("ctl./ctl."= dat[group_L=="control" & group_R=="control", log2FoldChange],
                "inact./inact."= dat[act_wilcox_L>=0.001 & act_wilcox_R>=0.001, log2FoldChange],
                "enh./ctl."= dat[act_wilcox_L<0.001 & group_R=="control", log2FoldChange],
                "ctl./enh."= dat[group_L=="control" & act_wilcox_R<0.001, log2FoldChange],
                "enh./enh."= dat[act_wilcox_L<0.001 & act_wilcox_R<0.001, log2FoldChange]), 
           violin= T, 
           compute_pval = list(c(1,2), c(1,3), c(2,4), c(3,4)),
           violcol = adjustcolor(c("lightgrey", "gold", "royalblue2", "royalblue2", "#74C27A"), 0.5),
           ylab= "Activity (log2)", 
           ylab.line = 2, 
           wilcox.alternative = "less")
abline(h= 0, lty= 2)
dev.off()
