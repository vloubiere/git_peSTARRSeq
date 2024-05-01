setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import and select active pairs
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
dat <- dat[(actPair)]

# Melt data for plotting ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat, 
           id.vars = c("L", "R", "Observed"), 
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]

# Check best predictor ----
t1 <- pl[, .(variable[which.min(abs(residuals))]), .(L, R)]$V1
t1 <- table(droplevels(t1))
names(t1) <- gsub(" model", "", names(t1))

# Plot ----
pdf("pdf/draft/Accuracy_predictions_pie_charts_with_lm.pdf",
    width = 3,
    height = 3)
par(mai= rep(1, 4), 
    mgp= c(1, .25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    lwd= .5)
pie(t1,
    labels = paste0(names(t1), " (", round(t1/sum(t1)*100), "%)\nn= ", formatC(t1, big.mark = ",")),
    cex= 8/12,
    lwd= .5,
    xpd= NA)
dev.off()