setwd("/groups/stark/vloubiere/projects/00006_vllib002_SCR1_validations_luc/")
source("/groups/stark/vloubiere/scripts/R_functions/my_plots.R")
require(data.table)

# Choose super additive and additive candidates
dat <- fread("data/F_exp_additive_scores_ctl_norm_0.25_no_DSCP_cutoff.txt")
dat[log2FC_obs>2, sa_cand_rank:= rank(-add_zscore)]
dat[log2FC_exp_L>1 & log2FC_exp_R>1, add_cand_rank:= rank(abs(add_zscore))]

# Chose negative controls
ctls <- fread("data/E_DESeq2_analysis_ctl_norm_0.25_no_DSCP_cutoff.txt")
ctls[grepl("control", enh_L), c("var_L", "median_L"):= .SD[grepl("control", enh_R), .(var(log2FC_obs, na.rm=T), median(log2FC_obs, na.rm=T))], enh_L]
ctls[grepl("control", enh_R), c("var_R", "median_R"):= .SD[grepl("control", enh_L), .(var(log2FC_obs, na.rm=T), median(log2FC_obs, na.rm=T))], enh_R]
L_ctl <- unique(na.omit(ctls[, .(enh_L, var_L, median_L)]))
R_ctl <- unique(na.omit(ctls[, .(enh_R, var_R, median_R)]))
ctls <- merge(L_ctl[, .(enh= enh_L, var_L, median_L)], R_ctl[, .(enh= enh_R, var_R, median_R)])
ctls[, var_rank:= rank(abs(median_L)+abs(median_R))]

# Assign (base on checks, see below)
ctls[var_rank %in% c(2,4,6,8), cdition:= "ctl_L"]
ctls[var_rank %in% c(5,7,9,10), cdition:= "ctl_R"]

# Checks
plot(ctls[, .(var_L, var_R)], col= ifelse(grepl("ctl", ctls$cdition), "red", "black"))
plot(ctls[, .(median_L, median_R)], col= ifelse(grepl("ctl", ctls$cdition), "red", "black"))
dds <- fread("data/E_DESeq2_analysis_ctl_norm_0.25_no_DSCP_cutoff.txt")
dds[ctls[var_rank<=16], neg_ctl_L:= i.var_rank, on= c("enh_L==enh")]
dds[ctls[var_rank<=16], neg_ctl_R:= i.var_rank, on= c("enh_R==enh")]
boxplot(log2FC_obs~neg_ctl_L, dds, notch= T)
boxplot(log2FC_obs~neg_ctl_R, dds, notch= T)

#### SUM UP TABLE AND SAVE
candidates <- data.table(dat[sa_cand_rank %in% 1:16, .(enh_L, enh_R, type= "super-additive")])
candidates <- rbind(candidates, data.table(dat[add_cand_rank %in% 1:16, .(enh_L, enh_R, type= "additive")]))
candidates <- rbind(candidates, cbind(ctls[cdition=="ctl_L", .(enh_L= enh)], ctls[cdition=="ctl_R", .(enh_R= enh, type= "ctl")]))
colnames(candidates) <- c("left", "right", "type")
candidates[, valid_ID:= paste0(substr(type, 1, 3), "_", .SD[,.I]), type]
candidates <- melt(candidates, id.vars = c("type", "valid_ID"))
colnames(candidates) <- c("cdition", "valid_ID", "position", "enh_ID")

fwrite(candidates, paste0("data/A_selected_candidates_", gsub("-", "", as.character(Sys.Date())), ".txt"))
