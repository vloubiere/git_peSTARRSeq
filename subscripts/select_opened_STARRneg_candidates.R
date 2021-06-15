dat <- data.table(file= c(list.files("../available_data_dm3/db/bed/", "ATAC", full.names = T),
                          list.files("../gw_STARRSeq_bernardo/db/raw_data/", "input_DSCP|DSCP.*Rep", full.names = T),
                          list.files("/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/mapped/", "mapped.reads.bed$", full.names = T)))

simp_expr <- "GSE119708_|_uniq.bed|_A.UMI_cut.bed|.UMI|_cut.bed|gw_|_gw_cut_merged.bed|S2_|.mapped.reads.bed"
dat[, name:= gsub(simp_expr, "", basename(file)), file]
dat[!grepl("_rep|_Rep", name), name:= paste0(name, "_rep1")]
dat[, variable:= gsub("rep|Rep", "", name)]