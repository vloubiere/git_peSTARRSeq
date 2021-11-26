setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import dat
dat <- list.files("db/final_tables_exp_model/counts_norm/", full.names = T)
names(dat) <- sapply(dat, function(x) paste0(unlist(tstrsplit(basename(x), "_", keep= c(1,3,4,5,6))), collapse= "_"))
dat <- rbindlist(lapply(dat, fread), idcol = "cdition", fill = T)
dat[grepl("_T8_", cdition), lib:= "ID_twist08"]
dat[grepl("_T12_", cdition), lib:= "ID_twist12"]
dat <- dat[!grepl("ts_", L) & !grepl("ts_", R)]

# Merge with features
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
feat <- melt(feat, measure.vars = c("ID_twist08", "ID_twist12"))[!is.na(value)]
dat <- merge(dat, 
             feat,
             by.x=c("L", "lib"),
             by.y=c("value", "variable"))
dat <- merge(dat, 
             feat,
             by.x=c("R", "lib"),
             by.y=c("value", "variable"),
             suffixes= c("_L", "_R"))
dat[, coor_L:= paste0(seqnames_L, ":", start_L, "-", end_L, ":", strand_L)]
dat[, coor_R:= paste0(seqnames_R, ":", start_R, "-", end_R, ":", strand_R)]
dat[, additive:= log2(2^median_L+2^median_R)]
dat[, diff:= log2FoldChange-additive]
dat[, additive_merge:= log2(2^median_merge_L+2^median_merge_R)]
dat[, diff_merge:= log2FC_merge-additive_merge]

# make groups as factors
dat[, group_L:= factor(group_L, 
                       levels = c("Silencer", "SUHW_peak", "control", "inducible", "OSC", "CP", "DHS_peak", "hk", "dev"))]
dat[, group_R:= factor(group_R, 
                       levels = c("Silencer", "SUHW_peak", "control", "inducible", "OSC", "CP", "DHS_peak", "hk", "dev"))]

# Add plot groups
# Twist 8 has only been used with DSCP
dat[grepl("_T8_", cdition)
    & L %in% dat[cdition=="vllib002_DSCP_T8_SCR1_300" & median_L>1, L]
    & group_L=="dev", active_group_L:= group_L]
dat[grepl("_T8_", cdition)
    & R %in% dat[cdition=="vllib002_DSCP_T8_SCR1_300" & median_R>1, R]
    & group_R=="dev", active_group_R:= group_R]
dat[grepl("_T8_", cdition)
    & active_group_L=="dev" 
    & active_group_R=="dev", active_group:= "dev"]
# Twist 12 has only been used with DSCP
dat[grepl("_DSCP_T12_", cdition)
    & group_L=="dev", active_group_L:= "dev"]
dat[grepl("_RpS12_T12_", cdition)
    & group_L=="hk", active_group_L:= "hk"]
dat[grepl("_DSCP_T12_", cdition)
    & group_R=="dev", active_group_R:= "dev"]
dat[grepl("_RpS12_T12_", cdition)
    & group_R=="hk", active_group_R:= "hk"]
dat[grepl("_T12_", cdition)
    & active_group_L %in% c("dev", "hk")
    & active_group_L==active_group_R, active_group:= active_group_L]

# CLean
cols <- c("cdition",
          "L", "R", 
          "group_L", "group_R", 
          "median_L", "median_R",
          "coor_L", "coor_R",
          "median_rep1_L", "median_rep2_L", "median_merge_L", "median_L",
          "median_rep1_R", "median_rep2_R", "median_merge_R", "median_R",
          "log2FC_rep1", "log2FC_rep2", "log2FC_merge", "log2FoldChange",
          "additive", "diff",
          "additive_merge", "diff_merge",
          "active_group_L", "active_group_R", "active_group",
          "dev_log2FoldChange_L", "dev_log2FoldChange_R",
          "hk_log2FoldChange_L", "hk_log2FoldChange_R",
          "col_L", "col_R",
          "closest_tss_L", "closest_tss_R",
          "closest_sel_tss_L", "closest_sel_tss_R",
          c(sort(grep("_log2FC_", colnames(dat), value = T))),
          c(sort(grep("^motif", colnames(dat), value = T))))
dat <- dat[, ..cols]

saveRDS(dat, "Rdata/final_results_table.rds")
