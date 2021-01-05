setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(yarrr)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
feat <- readRDS("Rdata/library/lib_features.rds")

#-------------------------------------------------#
# pdf("pdf/peSTARRSeq/smooth_scattreplots.pdf", 8, 9)
# par(mfrow= c(2,2), las= 1)
# act_cutoff <- 1
# 
# # Left enhancer individual act
# pl <- dat[median_R < act_cutoff, .(median_L, log2FoldChange)] 
# smoothScatter(pl, xlab= "X~ctl activity (left individual activity, log2)", ylab= "X~inactive activity (log2)")
# .lm <- lm(log2FoldChange~median_L, pl)
# abline(.lm, lty= 2)
# abline(0, 1)
# rsq <- round(summary(.lm)$r.square, 2)
# PCC <- round(cor(pl)[1,2], 2)
# legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)
# my_fig_label("A", cex= 2)
# 
# # Right enhancer individual act
# pl <- dat[median_L < act_cutoff, .(median_R, log2FoldChange)] 
# smoothScatter(pl, xlab= "ctl~X activity (right individual activity, log2)", ylab= "activity inactive~X (log2)")
# .lm <- lm(log2FoldChange~median_R, pl)
# abline(.lm, lty= 2)
# abline(0, 1)
# rsq <- round(summary(.lm)$r.square, 2)
# PCC <- round(cor(pl)[1,2], 2)
# legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)
# my_fig_label("B", cex= 2)
# 
# # Active enhancer pairs
# pl <- dat[median_L > act_cutoff & median_R > act_cutoff, .(log2FC_add, log2FoldChange)] 
# smoothScatter(pl, xlab= "additive score (log2)", ylab= "activity (log2)")
# .lm <- lm(log2FoldChange~log2FC_add, pl)
# abline(.lm, lty= 2)
# abline(0, 1)
# rsq <- round(summary(.lm)$r.square, 2)
# PCC <- round(cor(pl)[1,2], 2)
# legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)
# my_fig_label("C", cex= 2)
# 
# # Boxplot additive vs reality
# par(mar= c(5.1, 8.1, 4.1, 8.1))
# bpl <- melt(pl, measure.vars = colnames(pl))
# my_boxplot(value~variable, bpl, staplewex= 0, boxwex= 0.2, outline=F, notch= T, lty= 1, ylim= c(0, 10), 
#            ylab= "activity (log2)", names= c("add.", 'obs.'), pval_list = list(c(1,2)))
# my_fig_label("D", cex= 2)
# dev.off()

#----------------------------------------------#
pdf("pdf/peSTARRSeq/smooth_scatterplots_scaledxy.pdf", 8, 9)
par(mfrow= c(2,2), las= 1)
act_cutoff <- 1

ind_enh_lim = c(-3, 9);
paired_enh_lim = c(1, 9);

# Left enhancer individual act
pl <- dat[median_R < act_cutoff, .(median_L, log2FoldChange)] 
smoothScatter(pl, xlab= "X~ctl activity (left individual activity, log2)", ylab= "X~inactive activity (log2)", xlim= ind_enh_lim, ylim= ind_enh_lim)
.lm <- lm(log2FoldChange~median_L, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(pl)[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)
my_fig_label("A", cex= 2)

# Right enhancer individual act
pl <- dat[median_L < act_cutoff, .(median_R, log2FoldChange)] 
smoothScatter(pl, xlab= "ctl~X activity (right individual activity, log2)", ylab= "activity inactive~X (log2)", xlim= ind_enh_lim, ylim= ind_enh_lim)
.lm <- lm(log2FoldChange~median_R, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(pl)[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)
my_fig_label("B", cex= 2)

# Active enhancer pairs
pl <- dat[median_L > act_cutoff & median_R > act_cutoff, .(log2FC_add, log2FoldChange)] 
smoothScatter(pl, xlab= "additive score (log2)", ylab= "activity (log2)", xlim= paired_enh_lim, ylim= paired_enh_lim)
.lm <- lm(log2FoldChange~log2FC_add, pl)
abline(.lm, lty= 2)
abline(0, 1)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor(pl)[1,2], 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)
my_fig_label("C", cex= 2)

# Boxplot additive vs reality
par(mar= c(5.1, 8.1, 4.1, 8.1))
bpl <- melt(pl, measure.vars = colnames(pl))
my_boxplot(value~variable, bpl, staplewex= 0, boxwex= 0.2, outline=F, notch= T, lty= 1, ylim= paired_enh_lim, 
           ylab= "activity (log2)", names= c("add.", 'obs.'), pval_list = list(c(1,2)))
my_fig_label("D", cex= 2)
dev.off()
