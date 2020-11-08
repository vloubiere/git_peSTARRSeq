setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
lib <- readRDS("Rdata/library/lib_features.rds")
dat[lib, BA_L:= i.dev_log2FoldChange, on= "enh_L==ID"]
dat[lib, BA_R:= i.dev_log2FoldChange, on= "enh_R==ID"]
dat[, TWIST_additive_score:= log2(2^BA_L+2^BA_R)]


pdf("pdf/peSTARRSeq/BA_additive_scores_vs_peSTARRSeq.pdf")
smoothScatter(dat[, .(TWIST_additive_score, log2FoldChange)], las= 1)
abline(0, 1)
dev.off()