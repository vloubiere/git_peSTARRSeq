setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(yarrr)
require(ggplot2)
require(patchwork)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/all_expected_score.rds")
feat <- readRDS("Rdata/library/lib_features.rds")

#----------------------------------------------------------------#
# 1- format data
#----------------------------------------------------------------#
pl <- dat[!is.na(log2FoldChange) & !is.na(median_L) & !is.na(median_R)]

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
dmat <- dcast(pl, group_L~group_R, value.var = "log2FoldChange", fun.aggregate = median, na.rm= T)
dmat <- dmat[c(3,4,2,5,1)]
dmat <- dmat[, c(1,2,6,3,5,4)]

#-----------------------------------------------------#
# PLOT pheatmap
#-----------------------------------------------------#
pdf("pdf/peSTARRSeq/Heatmap_activity_type_aggregate.pdf")
layout(matrix(1:4, ncol= 2), widths= c(0.5,1), heights= c(0.46, 1))
Cc <- c("royalblue2", "gold", "lightgrey", "tomato", "#74C27A")
at <- c(3, 5, 4, 2, 1)
plot.new()

par(mar= c(4.5,5,1,1))
box <- my_boxplot(formula = -median_L~group_L, data= pl[, .(group_L, median_L), .(enh_L, group_L, median_L)], 
                  outline= T, col_box = Cc, horizontal= T, at= at, xaxt="n",
                  xlab= "left candidate\nindividual activity", xlim=c(0.7,5.3))
axis(1, at = seq(2, -6, -2), labels = seq(-2, 6, 2))
abline(v= box$DT_plot[.id=="control", lim3], lty= 2)

par(mar= c(1,1,4.5,5))
box <- my_boxplot(median_R~group_R, pl[, .(group_R, median_R), .(enh_R, group_R, median_R)], 
                  outline= T, col_box = Cc, at= at, xaxt= "n", xlim=c(0.7,5.3))
axis(3, at = c(3, 5, 4, 2, 1), labels = box$names)
mtext("right candidate\nindividual activity", side= 2, line = 2, cex= .8)
abline(h= box$DT_plot[.id=="control", lim3], lty= 2)

par(mar= c(4.5,1,1,5))
my_pheatmap(as.matrix(dmat, 1), cluster_rows = F, cluster_cols = F, row_labels = F, plot_grid = T, 
            display_numbers = T, legend_title = "median\nactivity\n(log2)", legend_adj_top = -0.1)
dev.off()


#-----------------------------------------------------#
# Violin plot categories
#-----------------------------------------------------#
sub <- pl[(group_L=="control" & group_R=="control") | group_L=="dev" | group_R=="dev"]
sub <- sub[enh_L!=enh_R & grepl("_C_", enh_L) & grepl("_C_", enh_R)]
# sub <- sub[grepl("_C_", enh_L) & grepl("_C_", enh_R)]

pdf("pdf/peSTARRSeq/Vioplot_activity_type_pairs_aggregate.pdf", height = 6)
par(mar= c(8,4,2,2))
pval <- list(c(1,2), c(2,3), c(2,4), c(2,5), c(6,7), c(6,8), c(6,9), c(6,10), c(2,6), c(2,10))
my_boxplot(formula = log2FoldChange~group_L+group_R, data = sub, las= 2, remove_empty = T, 
           main= "developmental enhancer containing pairs", 
           col_box = c("lightgrey", "#74C27A", rep("black", 3), "#74C27A",  rep("black", 3), "#74C27A"),
           at= c(1,2,6,10,9,8,7,5,4,3), pval_list = pval, ylim= c(-4, 10.5))
abline(h=0, lty=2)
dev.off()







