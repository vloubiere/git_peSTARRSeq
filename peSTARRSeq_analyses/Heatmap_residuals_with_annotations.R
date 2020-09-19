setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(yarrr)

dat <- readRDS("Rdata/SCR1_peSTARRSeq_final_table.rds")
feat <- readRDS("Rdata/lib_features.rds")

#----------------------------------------------------------------#
# Filtering
pl <- dat[!is.na(diff)]
pl[, check1:= .N>500, enh_L]
pl[, check2:= .N>500, enh_R]
pl <- pl[(check1) & (check2)] 

dmat <- dcast(pl, median_L+enh_L~median_R+enh_R, value.var = "diff")
dmat <- dmat[nrow(dmat):1]
mat <- 2^as.matrix(dmat[,-1], 1)

# Color codes
class_Cc <- data.table(class= c("hk", "dev", "OSC", "inducible"),
                       col= c("tomato", "#74C27A", "royalblue2", "#D3D3D3"),
                       value= 1:4)
setkeyv(class_Cc, "class")
# Cc <- colorRampPalette(c("cornflowerblue", "white", "white", "tomato"))
Cc <- colorRampPalette(c("blue", "yellow"))

# Plotting parameters
heat_lim <- c(0, 5)


#----------------------------------------------------------------#
pdf("pdf/Heatmap_residuals_with_annotation.pdf", 10, 10)
layout(matrix(1:16, ncol=4), widths= c(0.1,1,0.03,0.1), heights= c(0.1,1,0.03,0.1))

### ACTIVITY
# left
par(mar=c(0.25,0.25,0.25,0.25), xaxs="i", yaxs="i")
plot.new()
bar <- barplot(-rev(dmat$median_L), horiz = T, border= NA, col= ifelse(rev(dmat$median_L)<=0, "blue", "red"), axes=F, space = 0, xlim= c(-9,0.5))
segments(0, 0, 0, max(bar), lty= 1, col= "grey")
segments(-1, 0, -1, max(bar), lty= 2)
axis(3, at= c(-8, -1), labels = c(8, 1), line = 1)
# right
plot.new()
plot.new()
var <- as.numeric(unlist(tstrsplit(colnames(mat), "_", keep=1)))
bar <- barplot(var, border= NA, col= ifelse(var<=0, "blue", "red"), axes=F, space = 0, ylim= c(-0.5,9))
segments(0, 0, max(bar), 0, lty= 1, col= "grey")
segments(0, 1, max(bar), 1, lty= 2)
axis(2, at= c(1, 8), labels = c(1, 8), line = 1, las= 1)

### PAIR
my_pheatmap(mat, cluster_rows = F, cluster_cols = F, lim= heat_lim, col= Cc(100), plot_legend = F)

### RIGHT FEATURES
# Class
c_R <- t(as.matrix(class_Cc[unlist(tstrsplit(colnames(mat), "_", keep=2)), value]))
my_pheatmap(c_R, col= class_Cc$col, plot_grid = F, cluster_cols=F, plot_legend = F, col_labels = F, row_labels = F)
# Marginals
var <- log2(apply(mat, 2, median, na.rm= T))
barplot(-var, border= NA, col= ifelse(var<=0, "blue", "red"), axes=F, space = 0, ylim= c(-3, 3))
axis(4, at= 0:-2, labels = 0:2, line = 1, las= 1)

### LEFT FEATURES
# Class
plot.new()
c_L <- as.matrix(class_Cc[feat[match(dmat$enh_L, ID), group], value])
my_pheatmap(c_L, col= class_Cc$col, plot_grid = F, cluster_rows=F, plot_legend = F, col_labels = F, row_labels = F)
# Marginals
for(i in 1:3) plot.new()
var <- log2(apply(mat, 1, median, na.rm= T))
names(var) <- NULL
bar <- barplot(rev(var), horiz = T, border= NA, col= ifelse(rev(dmat$median_L)<=0, "blue", "red"), axes= F, space = 0, xlim= c(-2.5, 3))
segments(0, 0, 0, max(bar), lty= 1, col= "grey")
segments(1, 0, 1, max(bar), lty= 2)
axis(1, at= 0:2, labels = 0:2, line = 1, las= 1)

for(i in 1:7) plot.new()

# plot heatkey
par(mar= c(1,1,1,10))
my_pheatmap(mat, cluster_rows = F, cluster_cols = F, lim= heat_lim, col= Cc(100))

dev.off()
