setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table_high_cutoff.rds")
feat <- readRDS("Rdata/library/lib_features.rds")

# Left
.L <- unique(dat[, .(enh_L, median_L, sd_L)])
.L[feat, TWIST:= i.dev_log2FoldChange, on= "enh_L==ID"]
.R <- unique(dat[, .(enh_R, median_R, sd_R)])
.R[feat, TWIST:= i.dev_log2FoldChange, on= "enh_R==ID"]

.lmL <- lm(median_L~TWIST, .L)
.lmR <- lm(median_R~TWIST, .R)

pdf("pdf/peSTARRSeq/individual_scores_TWIST_compare_high_cutoff.pdf", width= 5, height = 9.5)
par(mfrow= c(2, 1), las= 1, mar= c(5,5,1,2))
col <- adjustcolor("grey", 0.5)
xl <- "TWIST STARR-Seq"

yl <- "X~ctl activity"
plot(.L[, .(TWIST, median_L)], xlab= xl, ylab= yl, pch= NA)
segments(.L$TWIST, .L$median_L-.L$sd_L, .L$TWIST, .L$median_L+.L$sd_L, col= col)
points(.L[, .(TWIST, median_L)], col= col, cex= 0.6, pch= 19)
points(.L[, .(TWIST, median_L)], cex= 0.6)
abline(.lmL)
rsq <- paste("R²", round(summary(.lmL)$r.squared, 2))
PCC <- paste("PCC=", round(cor.test(.L$median_L, .L$TWIST)$estimate, 2))
legend("topleft", legend = c(PCC, rsq), bty= "n")
my_fig_label("A", cex= 2)

yl <- "ctl~X activity"
plot(.R[, .(TWIST, median_R)], xlab= xl, ylab= yl, pch= NA)
segments(.R$TWIST, .R$median_R-.R$sd_R, .R$TWIST, .R$median_R+.R$sd_R, col= col)
points(.R[, .(TWIST, median_R)], col= col, cex= 0.6, pch= 19)
points(.R[, .(TWIST, median_R)], cex= 0.6)
abline(.lmR)
rsq <- paste("R²", round(summary(.lmR)$r.squared, 2))
PCC <- paste("PCC=", round(cor.test(.R$median_R, .R$TWIST)$estimate, 2))
legend("topleft", legend = c(PCC, rsq), bty= "n")
my_fig_label("B", cex= 2)
dev.off()
file.show("pdf/peSTARRSeq/individual_scores_TWIST_compare_high_cutoff.pdf")