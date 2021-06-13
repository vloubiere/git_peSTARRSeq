# Heatmap log2FC ####
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
dat <- dat[!(spike_in) & lib=="libvl013" & median_L>0.5 & median_R>0.5]
dat[, c("L", "R"):= .(substr(L, nchar(L)-4, nchar(L)), substr(R, nchar(R)-4, nchar(R)))]
dat[, diff:= log2FoldChange-add]
mat1 <- dcast(dat, -median_L+L~median_R+R, value.var = "log2FoldChange", fill = NA)
ind_L <- -(rev(mat1$median_L))
ind_R <- as.numeric(unlist(tstrsplit(colnames(mat1), "_", keep=1))[-c(1,2)])
colnames(mat1) <- unlist(tstrsplit(colnames(mat1), "_", keep = 2))
mat1 <- as.matrix(mat1[, -1], 1)
mat2 <- dcast(dat, -median_L+L~median_R+R, value.var = "diff", fill = NA)
mat2 <- as.matrix(mat2[, -c(1, 2)])
rownames(mat2) <- rownames(mat1)
colnames(mat2) <- colnames(mat1)

pdf("pdf/Heatmaps.pdf", 15, 15)
layout(matrix(1:4, ncol= 2), heights = c(0.25,1), widths = c(1,0.25))
#FC
par(mar= c(1,4.1,1,5.1), xaxs= "i", yaxs= "i")
bar <- barplot(ind_L)
points(bar, rev(apply(mat1, 1, mean, na.rm= T)))
par(mar= c(4.1,4.1,1,5.1))
my_heatmap(mat1, cluster_cols = F, cluster_rows = F, breaks = c(-2, 0, 10))
plot.new()
par(mar= c(4.1,1,1,1))
bar <- barplot(ind_R, horiz = T)
points(apply(mat1, 2, mean, na.rm= T), bar)
#diff
par(mar= c(1,4.1,1,5.1), xaxs= "i", yaxs= "i")
bar <- barplot(ind_L)
points(bar, rev(apply(mat2, 1, mean, na.rm= T)))
par(mar= c(4.1,4.1,1,5.1))
my_heatmap(mat2, cluster_cols = F, cluster_rows = F, breaks = c(-4, 0, 4))
plot.new()
par(mar= c(4.1,1,1,1))
bar <- barplot(ind_R, horiz = T)
points(apply(mat2, 2, mean, na.rm= T), bar)
dev.off()
####




