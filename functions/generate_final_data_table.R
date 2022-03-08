setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import dat
dat <- fread("Rdata/metadata_processed.txt")[(DESeq2)]
dat <- dat[, fread(FC_file), .(vllib, library, CP, spacer, spacer_size, FC_file)]
dat[, spacer:= paste0(spacer, "_", spacer_size), .(spacer, spacer_size)]
dat[, vllib:= factor(vllib, levels= sort(unique(dat$vllib)))]
dat[, library:= factor(library, c("T8", "T12"))]
dat$FC_file <- NULL
dat$spacer_size <- NULL
dat[, class_L:= ifelse(act_wilcox_L<0.05 & median_L>=1, "active", "inactive")]
dat[, class_R:= ifelse(act_wilcox_R<0.05 & median_R>=1, "active", "inactive")]

# Clean
saveRDS(dat, "Rdata/final_results_table.rds")
