setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import dat
dat <- fread("Rdata/metadata_processed.txt")[(DESeq2)]
dat <- dat[, fread(FC_file), .(vllib, CP, spacer, spacer_size, FC_file)]
dat[, spacer:= paste0(spacer, "_", spacer_size), .(spacer, spacer_size)]
dat$FC_file <- NULL
dat$spacer_size <- NULL

# Merge with features
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
feat[, coor:= paste0(seqnames, ":", start, "-", end, ":", strand)]
feat[, group:= factor(group, 
                      levels = c("repressor", "SUHW_peak", "control", "inducible", "OSC", "CP", "DHS_peak", "hk", "dev"))]
feat$oligo_full_sequence <- NULL
feat$linker_ID <- NULL
feat$seqnames <- NULL
feat$start <- NULL
feat$end <- NULL
feat$strand <- NULL
setkeyv(feat, "ID")
dat <- merge(dat,
             feat,
             by.x= "L",
             by.y= "ID")
dat <- merge(dat,
             feat,
             by.x= "R",
             by.y= "ID",
             suffixes= c("_L", "_R"))
dat[, additive:= log2(2^median_L+2^median_R)]
dat[, diff:= log2FoldChange-additive]

# Add Active column
dat[, active:= ifelse(median_L>1 & median_R>1, T, F)]

# Clean
setcolorder(dat,
            c("vllib", "CP", "L", "spacer", "R", 
              "group_L", "group_R", "detail_L", "detail_R",
              "median_L", "median_R", "log2FoldChange", "additive", "diff",
              "active", "col_L", "col_R", "coor_L", "coor_R",
              "motif_cl_L", "motif_cl_R",  
              "closest_tss_L", "closest_tss_R", "closest_sel_tss_L", "closest_sel_tss_R"))

saveRDS(dat, "Rdata/final_results_table.rds")
