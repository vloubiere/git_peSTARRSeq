setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Metadata ----
dat <- data.table(model= list.files("db/lasso_models/", "fold", full.names = T))
dat[, cdition:= gsub("_lasso_model.rds", "", basename(model))]
dat[, c("cdition", "var"):= .(gsub("(.*)_.*_.*$", "\\1", cdition),
                               gsub(".*_(.*$)", "\\1", cdition))]
dat[, rsq:= sapply(model, function(x) readRDS(x)$rsq)]
vl_par()
vl_boxplot(rsq~cdition+var,
           dat,
           tilt.names= T)

dat <- data.table(file= list.files("db/linear_models/", "^FC_lasso_.*predictions.rds$", full.names = T))
dat[, cdition:= tstrsplit(basename(file), "_", keep= 4)]
dat[, cdition:= switch(cdition, 
                       "ECD"= "S2 cells - ecdysone",
                       "large"= "S2 cells + ecdysone",
                       "OSC"= "OSC cells"), cdition]
dat <- dat[, readRDS(file)[, c("log2FoldChange", "residuals", "predicted_lasso", "predicted_residuals_lasso")], cdition]

# dcast before plotting ----
pl <- melt(dat,
            id.vars = "cdition",
            measure.vars = list(c("log2FoldChange", "residuals"),
                                c("predicted_lasso", "predicted_residuals_lasso")))
pl[, variable:= factor(c("Activity (log2)", "Residuals (log2)")[variable])]
setorderv(pl, c("cdition", "variable"))

# Plot ---
pdf("pdf/draft/residuals_prediction_lasso.pdf",
    height = 3, 
    width = 6)
par(mfrow= c(1,2),
    mai= rep(.9, 4), 
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
  smoothScatter(value1, 
                value2,
                xlab= "LASSO prediction (log2)",
                ylab= variable,
                col= adjustcolor(blues9[9], .3),
                main= cdition,
                xaxt= "n")
  axis(1, padj= -1.25)
  vl_plot_coeff(type = "rsq",
                value = vl_model_eval(value1, value2)$Rsquare,
                cex= 7/12,
                inset= c(-0.1, 0))
  .SD
}, .(cdition, variable)]
dev.off()