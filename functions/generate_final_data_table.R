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
dat[, class_L:= ifelse(act_wilcox_L<0.001, "active", "inactive")]
dat[, class_R:= ifelse(act_wilcox_R<0.001, "active", "inactive")]
dat[, class:= fcase(grepl("control", L), "ctl.", 
                    class_L=="active", "enh.",
                    class_L=="inactive", "inact.")]
dat[, class:= paste0(class, "/")]
dat[, class:= paste0(class, 
                     fcase(grepl("control", R), "ctl.", 
                           class_R=="active", "enh.",
                           class_R=="inactive", "inact."))]
dat[, class:= factor(class, c("ctl./ctl.",
                              "ctl./inact.",
                              "inact./ctl.",
                              "inact./inact.",
                              "ctl./enh.",
                              "inact./enh.",
                              "enh./ctl.",
                              "enh./inact.",
                              "enh./enh."))]
dat[, col:= c("grey0", "grey33", "grey66", "grey100", "cornflowerblue", "magenta1", "royalblue2", "darkorchid1", "#74C27A")[.GRP], keyby= class]
dat[, col:= adjustcolor(col, 0.5)]


# Clean
saveRDS(dat, "Rdata/final_results_table.rds")
