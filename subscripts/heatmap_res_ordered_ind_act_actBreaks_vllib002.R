setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")

# Matrix with activity breaks
mat <- dcast(dat,
             actL~actR,
             value.var = "residuals",
             fun.aggregate = mean)
mat <- as.matrix(mat, 1)
rownames(mat) <- dat[, paste0(actL, "\nn= ", length(unique(L))), actL]$V1
colnames(mat) <- dat[, paste0(actR, "\nn= ", length(unique(R))), actR]$V1
mat <- mat[nrow(mat):1,]
counts <- dcast(dat,
                actL~actR,
                value.var = "residuals",
                fun.aggregate = function(x) paste0("n= ", formatC(length(x), big.mark = ",")))
counts <- as.matrix(counts, 1)
counts <- counts[nrow(counts):1,]
hm <- vl_heatmap(mat, 
                 cluster_rows= F,
                 cluster_cols= F, 
                 breaks = c(-1, 0, 1),
                 legend_title= "Mean residuals (log2)",
                 plot= F,
                 tilt_colnames = T,
                 legend.cex = 0.6)

# plot ----
pdf("pdf/draft/heatmap_residuals_ordered_ind_act_actBreaks.pdf", 3.5, 2.4)
par(mar= c(3,4,1,5.5),
    mgp= c(2, 0.8, 0),
    tcl= -0.2,
    las= 1)
plot(hm)
dev.off()