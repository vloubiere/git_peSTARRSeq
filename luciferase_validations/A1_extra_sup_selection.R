setwd("/groups/stark/vloubiere/projects/0006_vllib002_SCR1_validations_luc/")
source("/groups/stark/vloubiere/scripts/R_functions/my_plots.R")
require(data.table)

dat <- readRDS("../pe_STARRSeq/Rdata/H4_additive_scores_merged_rep_high_cutoff.rds")
dat[, diff:= log2FoldChange-log2FC_add]

previous <- fread("Rdata/A_selected_candidates_20200521.txt")
dat <- dat[!enh_L %in% previous[position=="left", enh_ID] & !enh_R %in% previous[position=="right", enh_ID]]
dat <- dat[median_L> 2 & median_R>2]
dat <- dat[order(diff, decreasing= T)]

cand <- dat[c(1, 2, 7, 8, 13)]
cand[, valid_ID:= paste0("sup_", 17:21)]
cand[, cdition:= "super-additive"]
cand <- melt(cand, id.vars = c("cdition", "valid_ID"), measure.vars = c("enh_L", "enh_R"), variable.name = "position", value.name = "enh_ID")
cand[position=="enh_L", position:= "left"]
cand[position=="enh_R", position:= "right"]

fwrite(cand, "Rdata/A1_selected_extra_sup_20200625.txt")