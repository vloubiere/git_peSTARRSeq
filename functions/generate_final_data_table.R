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
# Add extra columns
dat[, additive:= log2(2^median_L+2^median_R)]
dat[, diff:= log2FoldChange-additive]
dat[, active:= ifelse(median_L>1 & median_R>1, T, F)]

# Add features
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- feat$add_feature(DT= dat, feature = feat$lib)

# Clean
saveRDS(dat, "Rdata/final_results_table.rds")
