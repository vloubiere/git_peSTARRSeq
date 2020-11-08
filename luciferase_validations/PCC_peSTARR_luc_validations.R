setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

dat <- readRDS("Rdata/luciferase_validations/C_luc_validations_final_table.rds")

#------------------------------------------------------------#
# 1- Add peSTARR-Seq data, groups and colors
#------------------------------------------------------------#
STARR <- readRDS("Rdata/processed_peSTARRSeq_data/DESeq2_FC_table.rds")
dat[STARR, log2FoldChange:= i.log2FoldChange, on= c("enh_L", "enh_R")]
feat <- readRDS("Rdata/library/lib_features.rds")
dat[, c("group_L", "col_L"):= feat[dat, .(group, col), on= "ID==enh_L"]]
dat[, c("group_R", "col_R"):= feat[dat, .(group, col), on= "ID==enh_R"]]
dat[, group := paste0(group_L, "~", group_R)]
dat[xor(group_L=="control", group_R=="control"), col:= ifelse(group_L=="control", col_R, col_L)]
dat[xor(group_L=="control", group_R=="control"), pch:= ifelse(group_L=="control", 17, 15)]
dat[!xor(group_L=="control", group_R=="control") & group_L!=group_R, col:= colorRampPalette(c(col_L, col_R))(3)[2], .(col_L, col_R)]
dat[!xor(group_L=="control", group_R=="control") & group_L==group_R, col:= col_L, .(col_L, col_R)]
dat[!xor(group_L=="control", group_R=="control"), pch:= 19]
dat <- dat[, .(luc_norm= mean(log2(luc_norm), na.rm= T), luc_sd= sd(log2(luc_norm), na.rm= T)), .(enh_L, enh_R, Sample_ID, log2FoldChange, group, col, pch)]

#------------------------------------------------------------#
# 2- PCC scatterplot
#------------------------------------------------------------#

pdf("pdf/luciferase_validations/PCC_luc_validations_peSTARRSeq.pdf", 6, 6.5)
par(las= 1)

plot(dat$log2FoldChange, dat$luc_norm,  col= dat$col, xlab= "pe-STARR-seq activity (log2)", ylab= "normalized luciferase activity (log2)", pch= dat$pch)
segments(dat$log2FoldChange, dat[, luc_norm-luc_sd], dat$log2FoldChange, dat[, luc_norm+luc_sd], col= dat$col)
.lm <- lm(luc_norm~log2FoldChange, dat)
abline(.lm, lty=2)
.lo <- loess(luc_norm~log2FoldChange, dat)
lines(seq(-5, 15, 0.1), predict(.lo, newdata = seq(-5, 15, 0.1)), lty=2, lwd=2, col= "red")
leg <- unique(dat[, .(group, col, pch)])
leg <- leg[order(col)]
legend("bottomright", legend = leg$group, col= leg$col, bty= "n", pch= leg$pch)
legend("topleft", legend = paste("R²=", round(summary(.lm)$r.squared, 2)), lty= 2, bty= "n")
dev.off()
