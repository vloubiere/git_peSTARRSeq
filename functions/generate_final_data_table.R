setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import dat
dat <- fread("Rdata/metadata_processed.txt")[(DESeq2)]
dat <- merge(unique(dat[, .(vllib, library, CP, spacer, spacer_size, FC_file)]),
             rbindlist(sapply(unique(dat$FC_file), fread), fill= T, idcol = "FC_file"),
             by= "FC_file")
dat <- dat[!is.na(log2FoldChange) & !is.na(median_L) & !is.na(median_R)]
dat[, spacer:= paste0(spacer, "_", spacer_size), .(spacer, spacer_size)]
dat[, vllib:= factor(vllib, levels= sort(unique(dat$vllib)))]
dat[, library:= factor(library, c("T8", "T12"))]
dat$FC_file <- NULL
dat$spacer_size <- NULL
dat[, class_L:= ifelse(act_wilcox_L<0.001 & median_L>log2(1.5), "active", "inactive")]
dat[is.na(class_L), class_L:= "inactive"]
dat[, class_R:= ifelse(act_wilcox_R<0.001 & median_R>log2(1.5), "active", "inactive")]
dat[is.na(class_R), class_R:= "inactive"]
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
                              "enh./ctl.",
                              "enh./inact.",
                              "ctl./enh.",
                              "inact./enh.",
                              "enh./enh."))]
dat[, col:= c("grey0", "grey33", "grey66", "grey100", "royalblue2", "cornflowerblue", "purple", "magenta", "#74C27A")[.GRP], keyby= class]
dat[, col:= adjustcolor(col, 0.5)]
test <- split(dat[vllib=="vllib015", log2FoldChange], 
              dat[vllib=="vllib015", .(col, class)], drop= T)
vl_boxplot(test, boxcol= unlist(tstrsplit(names(test), "[.]", keep= 1)), las= 2)

# Clean
saveRDS(dat, "Rdata/final_results_table.rds")
