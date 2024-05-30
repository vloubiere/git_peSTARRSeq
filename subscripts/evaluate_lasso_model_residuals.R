setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data ----
dat <- data.table(file= list.files("db/lasso_models/", "fold.*", full.names = T))
dat[, c("fold", "cdition"):= tstrsplit(basename(file), "_", keep= 1:2)]
dat[, cdition:= switch(cdition,
                       "large"= "S2 cells - ecdysone",
                       "ECD"= "S2 cells + ecdysone",
                       "OSC"= "OSC cells"), cdition]
dat <- dat[, .(model= .(readRDS(file))), .(cdition, fold)]

# Extract expected/observed ----
dat[, residuals:= lapply(model, function(x) x$data[set=="test", residuals])]
dat[, lasso_predicted_residuals:= lapply(model, function(x) x$data[set=="test", lasso_predicted_residuals])]
dat[, rsq:= sapply(model, function(x) x$residuals$adj.rsq)]
pl <- dat[, .(residuals= .(unlist(residuals)),
              lasso_predicted_residuals= .(unlist(lasso_predicted_residuals)),
              mean_rsq= paste0(round(mean(rsq), 2), "~", round(sd(rsq), 3))), cdition]

# Plot ---
pdf("pdf/draft/lasso_residuals_per_cdition.pdf",
    height = 3, 
    width = 9)
par(mfrow= c(1, 3),
    mai= rep(.9, 4), 
    mgp= c(1, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1,
    cex= 1)
pl[, {
  smoothScatter(lasso_predicted_residuals[[1]],
                residuals[[1]],
                xlab= "LASSO prediction (log2)",
                ylab= "Residuals (log2)",
                col= adjustcolor(blues9[9], .3),
                main= cdition,
                xaxt= "n")
  axis(1, padj= -1.25)
  vl_plot_coeff(type = "rsq",
                value = mean_rsq,
                cex= 7/12)
  .SD
}, cdition]
dev.off()