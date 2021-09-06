setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import features
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
.m <- melt(feat,
           measure.vars = list(c("ID_twist08", "ID_twist12"),
                               c("linker_twist08", "linker_twist12")), 
           value.name = c("ID", "linker"))
.m <- .m[, .(cdition= if(variable==1) 
  c("vllib002", "vllib006",  "vllib013", "vllib014") else if(variable==2)
    c("vllib015", "vllib016")), .m]
.m <- .m[!is.na(ID), !"variable"]

# Import dat
dat <- data.table(file= list.files("db/final_tables_exp_model/counts_norm/", full.names = T))
dat[, cdition:= tstrsplit(basename(file), "_", keep= 1)]
dat <- dat[, fread(file), (dat)]
dat <- dat[!grepl("ts_", L) & !grepl("ts_", R)]

# Merge
res <- merge(dat, 
             .m,
             by.x= c("cdition", "L"),
             by.y= c("cdition", "ID"),
             all.x= T)
res <- merge(res, 
             .m,
             by.x= c("cdition", "R"),
             by.y= c("cdition", "ID"),
             all.x= T, 
             suffixes= c("_L", "_R"))
res[, coor_L:= paste0(seqnames_L, ":", start_L, "-", end_L, ":", strand_L)]
res[, coor_R:= paste0(seqnames_R, ":", start_R, "-", end_R, ":", strand_R)]
res[, additive:= log2(2^median_L+2^median_R)]
res[, diff:= log2FoldChange-additive]

# Add plot groups
cmb <- CJ(res$group_L, res$group_R, unique = T)
cmb[, plot_group:= if(V1==V2)
  paste0(sort(c(V1, V2)), collapse = " x ") else
    paste0(paste0(sort(c(V1, V2)), collapse = " x "), " / ", paste0(rev(sort(c(V1, V2))), collapse = " x ")), .(V1, V2)]
cmb[, plot_group_LR:= paste0(V1, " x ", V2)]
res[cmb, c("plot_group", "plot_group_LR"):= .(i.plot_group, i.plot_group_LR), on= c("group_L==V1", "group_R==V2")]
res[cmb, c("plot_group", "plot_group_LR"):= .(i.plot_group, i.plot_group_LR), on= c("group_L==V2", "group_R==V1")]
res[, plot_group:= paste0(cdition, ": ", plot_group), .(cdition, plot_group)]
res[, plot_group_LR:= paste0(cdition, ": ", plot_group_LR), .(cdition, plot_group_LR)]
res[, active_plot_group:= as.character(NA)]
res[, active_plot_group_LR:= as.character(NA)]
res[cdition %in% c("vllib002", "vllib006", "vllib013", "vllib014", "vllib016") & 
      (group_L=="dev" | group_R=="dev"), c("active_plot_group", "active_plot_group_LR"):= .(plot_group, plot_group_LR)]
res[cdition == "vllib016" & 
      (group_L %in% c("CP", "dev", "hk") | 
         group_R %in% c("CP", "dev", "hk")), c("active_plot_group", "active_plot_group_LR"):= .(plot_group, plot_group_LR)]

# CLean
cols <- c("cdition",
          "L", "R", 
          "group_L", "group_R", 
          "median_L", "median_R",
          "coor_L", "coor_R", 
          "ctl_L", "ctl_R",
          "median_L", "median_R",
          "log2FC_rep1", "log2FC_rep2",
          "log2FoldChange",
          "additive",
          "diff",
          "plot_group", "active_plot_group",
          "plot_group_LR", "active_plot_group_LR",
          "active_plot_group",
          "dev_log2FoldChange_L", "dev_log2FoldChange_R",
          "hk_log2FoldChange_L", "hk_log2FoldChange_R",
          "col_L", "col_R",
          "closest_tss_L", "closest_tss_R",
          "closest_sel_tss_L", "closest_sel_tss_R",
          c(sort(grep("_log2FC_", colnames(res), value = T))),
          c(sort(grep("^motif__", colnames(res), value = T))))
res <- res[, ..cols]

saveRDS(res, "Rdata/final_results_table.rds")
