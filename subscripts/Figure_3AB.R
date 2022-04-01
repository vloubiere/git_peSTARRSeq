setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
lib <- readRDS("Rdata/final_results_table.rds")
lib[, diff:= log2FoldChange-additive]
lib <- feat$add_feature(lib, feat$lib)
dat <- merge(lib[vllib=="vllib015" & group_L %in% c("hk", "dev") & group_R %in% c("hk", "dev"), !c("vllib", "CP", "spacer")],
             lib[vllib=="vllib016" & group_L %in% c("hk", "dev") & group_R %in% c("hk", "dev"), !c("vllib", "CP", "spacer")],
             by= c("L", "R", "group_L", "group_R", "class_L", "class_R", "col"), 
             suffixes= c("_dev", "_hk"))
dat[group_L=="dev" & group_R=="dev", col:= "#74C27A"]
dat[group_L=="dev" & group_R=="hk", col:= "cyan"]
dat[group_L=="hk" & group_R=="dev", col:= "royalblue2"]
dat[group_L=="hk" & group_R=="hk", col:= "tomato"]

.m <- melt(dat, 
           id.vars = c("group_L", "group_R", "col"), 
           measure.vars = list("log2FoldChange_dev", "log2FoldChange_hk"))
dev_act <- split(.m[, .(value1, col)], .m[, .(group_L, group_R)])
dev_act <- dev_act[c("dev.dev", "dev.hk", "hk.dev", "hk.hk")]
hk_act <- split(.m[, .(value2, col)], .m[, .(group_L, group_R)])
hk_act <- hk_act[c("hk.hk", "hk.dev", "dev.hk", "dev.dev")]

.m <- melt(dat, 
           id.vars = c("group_L", "group_R", "col"), 
           measure.vars = list("diff_dev", "diff_hk"))
dev_diff <- split(.m[, .(value1, col)], .m[, .(group_L, group_R)])
dev_diff <- dev_diff[c("dev.dev", "dev.hk", "hk.dev", "hk.hk")]
hk_diff <- split(.m[, .(value2, col)], .m[, .(group_L, group_R)])
hk_diff <- hk_diff[c("hk.hk", "hk.dev", "dev.hk", "dev.dev")]

pdf("pdf/draft/Figure_3AB.pdf", 
    height = 5.5,
    width= 5.6)
par(las= 1,
    mar= c(5.1,4.1,0,0))
layout(matrix(c(2,4,1,3), ncol = 2, byrow = T), 
       widths = c(2,0.7),
       heights = c(0.62,2))
yl <- c(-2,8.5)
xl <- c(-2,12)
plot(dat$log2FoldChange_dev,
     dat$log2FoldChange_hk, 
     xlab= "Activity dCP (log2)",
     ylab= "Activity hkCP (log2)",
     pch= 19,
     col= adjustcolor(dat$col, 0.5),
     cex= 0.5,
     xlim= xl,
     ylim= yl)
abline(h= 0, lty= 2)
abline(v= 0, lty= 2)
abline(0,1,lty= 2)

par(mar= c(1,4.1,0.5,0))
vl_boxplot(lapply(dev_act, `[[`, 1), 
           ylim= xl,
           horizontal= T,
           xaxt= "n",
           violcol= adjustcolor(unlist(sapply(dev_act, function(x) x[1,2])), 0.5), 
           compute_pval = list(c(1,4)), 
           violin= T)
abline(v= 0, lty= 2)

par(mar= c(5.1,1,0,0.5),
    las= 2)
vl_boxplot(lapply(hk_act, `[[`, 1), 
           ylim= yl,
           yaxt= "n",
           violcol= adjustcolor(unlist(sapply(hk_act, function(x) x[1,2])), 0.5), 
           compute_pval = list(c(1,4)),
           violin= T)
abline(h= 0, lty= 2)

#---------#
plot.new()
par(las= 1,
    mar= c(5.1,4.1,0,0))
yl <- c(-5,4)
xl <- c(-4,5)
plot(dat$diff_dev,
     dat$diff_hk, 
     xlab= "Observed/Exp. Add. dCP (log2)",
     ylab= "Observed/Exp. Add. hkCP (log2)",
     pch= 19,
     col= adjustcolor(dat$col, 0.5),
     cex= 0.5,
     xlim= xl,
     ylim= yl)
abline(h= 0, lty= 2)
abline(v= 0, lty= 2)
abline(0,1,lty= 2)

par(mar= c(1,4.1,0.5,0))
vl_boxplot(lapply(dev_diff, `[[`, 1), 
           ylim= xl,
           horizontal= T,
           xaxt= "n",
           violcol= adjustcolor(unlist(sapply(dev_diff, function(x) x[1,2])), 0.5), 
           compute_pval = list(c(1,4)), 
           violin= T)
abline(v= 0, lty= 2)

par(mar= c(5.1,1,0,0.5),
    las= 2)
vl_boxplot(lapply(hk_diff, `[[`, 1), 
           ylim= yl,
           yaxt= "n",
           violcol= adjustcolor(unlist(sapply(hk_diff, function(x) x[1,2])), 0.5), 
           compute_pval = list(c(1,4)),
           violin= T)
abline(h= 0, lty= 2)
dev.off()




