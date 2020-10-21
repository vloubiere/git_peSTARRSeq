setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions", ".R$", full.names = T), source)
require(data.table)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/all_expected_score.rds")
feat <- readRDS("Rdata/library/lib_features.rds")

#----------------------------------------------------------------#
# 1- format data
#----------------------------------------------------------------#
pl <- dat[!is.na(log2FoldChange) & !is.na(median_L) & !is.na(median_R)]
pl[, check1:= .N>0, enh_L]
pl[, check2:= .N>0, enh_R]
pl <- pl[(check1) & (check2)]

# Add group
pl[feat, group_L:= i.group, on= "enh_L==ID"]
pl[feat, group_R:= i.group, on= "enh_R==ID"]
# Handle E coli cases
pl[is.na(group_L), group_L:="control"]
pl[is.na(group_R), group_R:="control"]

# Add group color
class_Cc <- data.table(class= c("hk", "dev", "OSC", "inducible", "control"),
                       col= c("tomato", "#74C27A", "royalblue2", "gold", "#D3D3D3"))
pl[class_Cc, col_L:= i.col, on= "group_L==class"]
pl[class_Cc, col_R:= i.col, on= "group_R==class"]

# dcast/melt
pl[, Rvar:= paste0("_", paste(c(enh_R, group_R, col_R), collapse = "__")), .(enh_R, group_R, col_R)]
dmat <- dcast(pl, -median_L+enh_L+group_L+col_L~median_R+Rvar, value.var = "log2FoldChange", fill = NA)
mat <- as.matrix(dmat[, -c(1:4)])
pl <- melt(dmat, id.vars = c("median_L", "enh_L", "group_L", "col_L"), value.name = "log2FoldChange", variable.name = "Rvar")
pl[, median_L:= -median_L]
pl[, c("median_R", "enh_R", "group_R", "col_R"):= tstrsplit(Rvar, "__")]
pl$Rvar <- NULL

#----------------------------------------------------------------#
# 2- plot
#----------------------------------------------------------------#
# pdf
pdf("pdf/peSTARRSeq/Heatmap_activity_with_annotation.pdf", 13, 11.5)
layout(matrix(1:16, ncol=4), widths= c(0.1,0.03,0.1,1), heights= c(0.1,0.03,0.1,1))
par(mar=c(0.25,0.5,0.5,0.25), xaxs="i", yaxs="i")

####-------LEFT INDIVIDUAL ACTIVITY, FEATURES and MARGINALS-------####

# median
.c <- as.numeric(pl[, median_L, .(enh_L, median_L)]$median_L)
plot.new()
plot.new()
plot.new()
bar <- barplot(-rev(.c), horiz = T, border= NA, col= ifelse(rev(.c)<=0, "blue", "red"), axes=F, space = 0, xlim= c(-9,0.5))
segments(0, 0, 0, max(bar), lty= 1, col= "grey")
segments(-1, 0, -1, max(bar), lty= 2)
axis(3, at= c(-8, -1), labels = c(8, 1), line = 1)

# group
.c <- rev(pl[, col_L, enh_L]$col_L)
plot.new()
plot.new()
plot.new()
barplot(rep(1, length(.c)), horiz = T, border= NA, col= .c, axes=F, space = 0)

# Marginals
.c <- rev(pl[, median(log2FoldChange, na.rm= T), enh_L]$V1)
plot.new()
plot.new()
plot.new()
bar <- barplot(-.c, horiz = T, border= NA, col= ifelse(.c<=0, "blue", "red"), axes= F, space = 0, xlim= c(-7, 0))
segments(0, 0, 0, max(bar), lty= 1, col= "grey")
segments(-1, 0, -1, max(bar), lty= 2)
axis(3, at= c(-1, -7), labels = c(1, 7), line = 1, las= 1)

####-------RIGHT INDIVIDUAL ACTIVITY, FEATURES and MARGINALS-------####
par(mar= c(0.25,0.5,0.5,12), xaxs="i", yaxs="i")
# median
.c <- as.numeric(pl[, median_R, .(enh_R, median_R)]$median_R)
bar <- barplot(.c, border= NA, col= ifelse(.c <= 0, "blue", "red"), axes=F, space = 0, ylim= c(-0.5,9))
segments(0, 0, max(bar), 0, lty= 1, col= "grey")
segments(0, 1, max(bar), 1, lty= 2)
axis(2, at= c(1, 8), labels = c(1, 8), line = 1, las= 1)

# group
.c <- pl[, col_R, enh_R]$col_R
barplot(rep(1, length(.c)), border= NA, col= .c, axes=F, space = 0)

# Marginals
.c <- pl[, median(log2FoldChange, na.rm= T), enh_R]$V1
bar <- barplot(.c, border= NA, col= ifelse(.c<=0, "blue", "red"), axes=F, space = 0, ylim= c(0, 7))
segments(0, 0, max(bar), 0, lty= 1, col= "grey")
segments(0, 1, max(bar), 1, lty= 2)
axis(2, at= c(1, 7), labels = c(1, 7), line = 1, las= 1)

####-------ADDITIVE SCORES-------####
# PAIR
my_pheatmap(mat, cluster_rows = F, cluster_cols = F, breaks= c(-2,0,9),
            col= c("blue", "white", "red"), bg_col= "black", legend_cex = 1.8,
            legend_title= "Activity\n(log2)", legend_adj_top = -0.025, legend_adj_left = -0.025)

####-------LEGEND-------####
points(1.055, 0.675, pch=15, xpd= T, cex= 3, col= "black")
text(1.055, 0.675, pos= 4, "NA", xpd= T, offset = 1, cex= 1.5)
points(rep(1.05, 5), seq(0.425, 0.55, length.out = 5), pch=15, xpd= T, cex= 3, col= class_Cc$col)
text(rep(1.05, 5), seq(0.425, 0.55, length.out = 5), pos= 4, class_Cc$class, xpd= T, offset = 1, cex= 1.5)
text(1.033, 0.6, "Enhancer\ntype", xpd= T, cex= 1.8, pos= 4)

dev.off()
