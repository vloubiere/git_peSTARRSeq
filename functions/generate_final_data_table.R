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
.m[group=="shared", group:= ifelse(dev_log2FoldChange>hk_log2FoldChange, "dev", "hk")]

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
cols <- c("cdition", 
          "L", "R", 
          "group_L", "group_R", 
          "median_L", "median_R",
          "coor_L", "coor_R", 
          "ctl_L", "ctl_R",
          "median_L", "median_R",
          "log2FoldChange",
          "additive",
          "diff",
          "dev_log2FoldChange_L", "dev_log2FoldChange_R",
          "hk_log2FoldChange_L", "hk_log2FoldChange_R",
          "col_L", "col_R",
          "closest_tss_L", "closest_tss_R",
          "closest_sel_tss_L", "closest_sel_tss_R",
          c(sort(grep("_log2FC_", colnames(res), value = T))),
          c(sort(grep("^motif__", colnames(res), value = T))))
res <- res[, ..cols]


saveRDS(res, "Rdata/final_results_table.rds")
