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
dat[, group:= paste0(feat[dat, group, on= "ID==enh_L"], "~", feat[dat, group, on= "ID==enh_R"])]
dat[, group:= gsub("NA", "control", group)]
Cc <- data.table(group= c("hk~hk", "control~dev", "control~hk", "dev~control", "control~control", "dev~dev", "hk~dev", "hk~control", "dev~hk"),
                 col= c("red", "cornflowerblue", "chocolate1", "cyan3", "grey", "darkblue", "gold", "coral4", "green"))
dat[Cc, col:= i.col, on= "group"]

#------------------------------------------------------------#
# 2- PCC scatterplot
#------------------------------------------------------------#

pdf("pdf/luciferase_validations/PCC_luc_validations_peSTARRSeq.pdf", 6, 6.5)
par(las= 1)
pl <- dat[, .(log2FoldChange[1], mean(log2(luc_norm)), 
              mean(log2(luc_norm))+sd(log2(luc_norm)), mean(log2(luc_norm))-sd(log2(luc_norm)),  col), Sample_ID]

plot(pl$V1, pl$V2, col= pl$col, pch= 19, xlab= "pe-STARR-seq activity (log2)", ylab= "normalized luciferase activity (log2)")
segments(pl$V1, pl$V3, pl$V1, pl$V4, col= pl$col)
.lm <- lm(V2~V1, pl)
abline(.lm, lty=2)
leg <- unique(dat[, .(group, col)])
legend("bottomright", legend = leg$group, col= leg$col, bty= "n", pch= 19)
legend("topleft", legend = paste("R²=", round(summary(.lm)$r.squared, 2)), lty= 2, bty= "n")
dev.off()









