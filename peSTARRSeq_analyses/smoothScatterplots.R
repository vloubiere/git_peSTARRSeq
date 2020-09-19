setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(yarrr)

dat <- readRDS("Rdata/SCR1_peSTARRSeq_final_table.rds")
feat <- readRDS("Rdata/lib_features.rds")

#-------------------------------------------------#
pdf("pdf/smooth_scattreplots.pdf", 16, 4)
layout(matrix(1:5, ncol=5), widths = c(1, 1, 1, 1, 0.5))
par(las= 1)
act_cutoff <- 1

# Direction effect
pl <- dat[, .(log2FoldChange, log2FoldChange_rev)] 
smoothScatter(pl, xlab= "X~Y activity (log2)", ylab= "Y~X activity (log2)")
.lm <- lm(log2FoldChange_rev~log2FoldChange, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(na.omit(pl))[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)

# Left enhancer individual act
pl <- dat[median_R < act_cutoff, .(median_L, log2FoldChange)] 
smoothScatter(pl, xlab= "median X~ctl (log2)", ylab= "activity X~inactive (log2)")
.lm <- lm(log2FoldChange~median_L, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(pl)[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)

# Right enhancer individual act
pl <- dat[median_L < act_cutoff, .(median_R, log2FoldChange)] 
smoothScatter(pl, xlab= "median ctl~X (log2)", ylab= "activity inactive~X (log2)")
.lm <- lm(log2FoldChange~median_R, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(pl)[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)

# Active enhancer pairs
pl <- dat[median_L > act_cutoff & median_R > act_cutoff, .(log2FC_add, log2FoldChange)] 
smoothScatter(pl, xlab= "additive score (log2)", ylab= "activity (log2)")
.lm <- lm(log2FoldChange~log2FC_add, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(pl)[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)

# Boxplot additive vs reality
bpl <- melt(pl, measure.vars = colnames(pl))
boxplot(value~variable, bpl, staplewex= 0, boxwex= 0.2, outline=F, notch= T, lty= 1, ylim= c(0, 11.5), 
        ylab= "activity (log2)", names= c("add.", 'obs.'))
my_vioplot(bpl, "value", "variable")
my_pval_plot(value~variable, bpl, y_quant = 0.9999)
dev.off()
