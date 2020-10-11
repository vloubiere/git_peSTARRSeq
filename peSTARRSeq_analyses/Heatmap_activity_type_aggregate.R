setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(yarrr)
require(ggplot2)
require(patchwork)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/expected_score.rds")
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
box <- my_boxplot(-median_L~group_L, pl[, .(group_L, median_L), .(enh_L, group_L, median_L)], 
                  outline= T, col_box = Cc, horizontal= T, at= at, xaxt="n",
                  xlab= "left candidate\nindividual activity", xlim=c(0.7,5.3))
axis(1, at = seq(2, -6, -2), labels = seq(-2, 6, 2))
abline(v= box$stats[3,3], lty= 2)

par(mar= c(1,1,4.5,5))
box <- my_boxplot(median_R~group_R, pl[, .(group_R, median_R), .(enh_R, group_R, median_R)], 
                  outline= T, col_box = Cc, at= at, xaxt= "n", xlim=c(0.7,5.3))
axis(3, at = c(3, 5, 4, 2, 1), labels = box$names)
mtext("right candidate\nindividual activity", side= 2, line = 2, cex= .8)
abline(h= box$stats[3,3], lty= 2)

par(mar= c(4.5,1,1,5))
my_pheatmap(as.matrix(dmat, 1), cluster_rows = F, cluster_cols = F, row_labels = F, plot_grid = T, 
            display_numbers = T, legend_title = "median\nactivity\n(log2)", legend_adj_top = -0.1)
dev.off()


#-----------------------------------------------------#
# Violin plot categories
#-----------------------------------------------------#

pdf("pdf/peSTARRSeq/Vioplot_activity_type_pairs_aggregate.pdf", width = 10)
par(mar= c(10,4,2,2))
at <- c(1,5,4,3,2,25,29,28,26,27,19,23,22,20,21,13,17,16,14,15,7,11,10,8,9)
Cc1 <- c("lightgrey", "royalblue2", "gold", "tomato", "#74C27A")
Cc2 <- c(sapply(Cc1, function(x) c(x, rep("black", 4))))
Cc2[1:5] <- Cc1
my_boxplot(formula = log2FoldChange~group_L+group_R, data = pl, las= 2, 
           at= at, ylim= c(-4, 11), col_box = Cc2,
           pval_list = list(c(1,2), c(1,3), c(1,4), c(1,5),
                            c(6,7), c(6,8), c(6,9), c(6,10),
                            c(11,12), c(11,13), c(11,14), c(11,15),
                            c(16,17), c(16,18), c(16,19), c(16,20),
                            c(21,22), c(21,23), c(21,24), c(21,25)))
abline(h=0, lty=2)
legend("topleft", fill= Cc1, legend= c("control", "OSC", "inducible", "hk", "dev"), bty= "n")
dev.off()







