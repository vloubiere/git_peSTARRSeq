setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------------#
# Import dat
#-----------------------------------------------------#
dat <- fread("Rdata/metadata_processed.txt")[(DESeq2)]

#-----------------------------------------------------#
# Clean
#-----------------------------------------------------#
dat <- merge(unique(dat[, .(vllib, library, CP, spacer, spacer_size, FC_file)]),
             rbindlist(sapply(unique(dat$FC_file), fread), fill= T, idcol = "FC_file"),
             by= "FC_file")
dat <- dat[!is.na(log2FoldChange) & !is.na(median_L) & !is.na(median_R)]
dat[, spacer:= paste0(spacer, "_", spacer_size), .(spacer, spacer_size)]
dat[, vllib:= factor(vllib, levels= sort(unique(vllib)))]
dat[, library:= factor(library, c("T8", "T12"))]
dat$FC_file <- dat$spacer_size <- NULL

#-----------------------------------------------------#
# Define active/inactive individual enhancers
#-----------------------------------------------------#
dat[, class_act_L:= fcase(FDR_L<0.05 & median_L>1, "active", default= "inactive")]
dat[, class_act_R:= fcase(FDR_R<0.05 & median_R>1, "active", default= "inactive")]

#-----------------------------------------------------#
# Activity classes -> cannot be defined without screen data in 300bp_uniq_enhancers
#-----------------------------------------------------#
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
# Activity classes colors
dat[, col_act:= c("grey0", 
                  "grey33", 
                  "grey66", 
                  "grey100", 
                  "royalblue2", 
                  "cornflowerblue", 
                  "purple", 
                  "magenta", 
                  "#74C27A")[.GRP], keyby= class_act]

# Clean
saveRDS(dat,
        "Rdata/final_results_table.rds")
