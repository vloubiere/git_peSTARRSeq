setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(vioplot)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
lib[, diff:= log2FoldChange-multiplicative]
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
# dat <- lib[vllib=="vllib002" & act_wilcox_L<0.001 & act_wilcox_R<0.001]
dat <- lib[vllib=="vllib002"]
dat <- feat$add_feature(dat, feat$lib)
dat <- feat$add_feature(dat, feat$closest_TSSs)


enh <- unique(dat[, c(L, R)])
pairs <- data.table(A= enh)
pairs <- pairs[, .(B= enh[-(1:.GRP)]), A]
pairs[dat, log2FoldChange.x:= i.log2FoldChange, on= c("A==L", "B==R")]
pairs[dat, log2FoldChange.y:= i.log2FoldChange, on= c("A==R", "B==L")]
pairs[dat, diff.x:= i.diff, on= c("A==L", "B==R")]
pairs[dat, diff.y:= i.diff, on= c("A==R", "B==L")]
pairs <- na.omit(pairs)


smoothScatter(pairs[, .(log2FoldChange.x, log2FoldChange.y)])
smoothScatter(pairs[, .(diff.x, diff.y)])
pairs[identify(pairs[, .(log2FoldChange.x, log2FoldChange.y)])]
points(pairs[grepl("dev", A) & grepl("dev", B), .(diff.x, diff.y)])
points(pairs[grepl("hk", A) | grepl("hk", B), .(diff.x, diff.y)], col= "red")


smoothScatter(dat[act_wilcox_L<0.05 & median_L>=1 & act_wilcox_R<0.05 & median_R>=1, .(additive, log2FoldChange)])
cor.test(dat[act_wilcox_L<0.05 & median_L>=1 & act_wilcox_R<0.05 & median_R>=1, additive],
         dat[act_wilcox_L<0.05 & median_L>=1 & act_wilcox_R<0.05 & median_R>=1, log2FoldChange])
abline(0,1)
