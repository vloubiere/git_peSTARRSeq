setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data ----
dat <- data.table(file= list.files("db/lasso_models/", "fold.*large", full.names = T))
dat[, fold:= tstrsplit(basename(file), "_", keep= 1)]
dat <- dat[, .(model= .(readRDS(file))), fold]

# Extract expected/observed ----
dat[, log2FoldChange:= lapply(model, function(x) x$data[set=="test", log2FoldChange])]
dat[, lasso_predicted_log2FoldChange:= lapply(model, function(x) x$data[set=="test", lasso_predicted_log2FoldChange])]
dat[, rsq:= sapply(model, function(x) x$log2FoldChange$adj.rsq)]
pl <- dat[, .(log2FoldChange= .(unlist(log2FoldChange)),
              lasso_predicted_log2FoldChange= .(unlist(lasso_predicted_log2FoldChange)),
              mean_rsq= paste0(round(mean(rsq), 2), "~", round(sd(rsq), 2)))]

# Plot ---
pdf("pdf/draft/lasso_large_WT_activity_prediction.pdf",
    height = 3, 
    width = 3)
par(mai= rep(.9, 4), 
    mgp= c(1, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
pl[, {
  smoothScatter(lasso_predicted_log2FoldChange[[1]],
                log2FoldChange[[1]],
                xlab= "LASSO prediction (log2)",
                ylab= "Activity (log2)",
                col= adjustcolor(blues9[9], .3),
                xaxt= "n")
  axis(1, padj= -1.25)
  vl_plot_coeff(type = "rsq",
                value = mean_rsq,
                adjusted = T,
                cex= 7/12,
                inset= c(-0.1, 0))
}]
dev.off()