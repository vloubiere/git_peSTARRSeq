setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import dat
dat <- fread("Rdata/metadata_processed.txt")[(DESeq2)]
# Clean
dat <- merge(unique(dat[, .(vllib, library, CP, spacer, spacer_size, FC_file)]),
             rbindlist(sapply(unique(dat$FC_file), fread), fill= T, idcol = "FC_file"),
             by= "FC_file")
dat <- dat[!is.na(log2FoldChange) & !is.na(median_L) & !is.na(median_R)]
dat[, spacer:= paste0(spacer, "_", spacer_size), .(spacer, spacer_size)]
dat[, vllib:= factor(vllib, levels= sort(unique(dat$vllib)))]
dat[, library:= factor(library, c("T8", "T12"))]
dat$FC_file <- NULL
dat$spacer_size <- NULL
# Define active/inactive individual enhancers
dat[, class_act_L:= ifelse(act_wilcox_L<0.001 & median_L>log2(1.5), "active", "inactive")]
dat[is.na(class_act_L), class_act_L:= "inactive"]
dat[, class_act_R:= ifelse(act_wilcox_R<0.001 & median_R>log2(1.5), "active", "inactive")]
dat[is.na(class_act_R), class_act_R:= "inactive"]
# Activity classes and colors
dat[, class_act:= fcase(grepl("control", L), "ctl.", 
                        class_act_L=="active", "enh.",
                        class_act_L=="inactive", "inact.")]
dat[, class_act:= paste0(class_act, "/")]
dat[, class_act:= paste0(class_act, 
                         fcase(grepl("control", R), "ctl.", 
                               class_act_R=="active", "enh.",
                               class_act_R=="inactive", "inact."))]
dat[, class_act:= factor(class_act, 
                         c("ctl./ctl.",
                           "ctl./inact.",
                           "inact./ctl.",
                           "inact./inact.",
                           "enh./ctl.",
                           "enh./inact.",
                           "ctl./enh.",
                           "inact./enh.",
                           "enh./enh."))]
dat[, col_act:= c("grey0", "grey33", "grey66", "grey100", "royalblue2", "cornflowerblue", "purple", "magenta", "#74C27A")[.GRP], keyby= class_act]
dat[, col_act:= c("grey0", "grey33", "grey66", "grey100", "royalblue2", "cornflowerblue", "purple", "magenta", "#74C27A")[class_act]]
test <- split(dat[vllib=="vllib015", log2FoldChange], 
              dat[vllib=="vllib015", .(col_act, class_act)], drop= T)
vl_boxplot(test, boxcol= unlist(tstrsplit(names(test), "[.]", keep= 1)), las= 2)
# hk/dev classes and colors
dat[grepl("dev", L) & grepl("dev", R), c("col_hkdev", "class_hkdev"):= .("#74C27A", "dev/dev")]
dat[grepl("dev", L) & grepl("hk", R), c("col_hkdev", "class_hkdev"):= .("cyan", "dev/hk")]
dat[grepl("hk", L) & grepl("dev", R), c("col_hkdev", "class_hkdev"):= .("royalblue2", "hk/dev")]
dat[grepl("hk", L) & grepl("hk", R), c("col_hkdev", "class_hkdev"):= .("tomato", "hk/hk")]

# Clean
saveRDS(dat, "Rdata/final_results_table.rds")