setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(readxl)

#------------------------------------------------------------#
# 1- Import and format data
#------------------------------------------------------------#
# Plate schemes
scheme <- data.table(file= list.files("db/luciferase/peSTARRSeq_validations/", "scheme", recursive = T, full.names = T))
scheme[, c("date", "plate") := tstrsplit(basename(file), "_|[.]", keep= c(1,3)), file]
scheme <- scheme[, melt(fread(file, 
                              header = T, 
                              colClasses = "character", 
                              na.strings = "", 
                              fill= T), 
                        id.vars = "rownames", 
                        value.name = "Sample_ID"), (scheme)]
colnames(scheme)[4:5] <- c("row", "col")
# Add technical replicates
rep <- fread("db/luciferase/peSTARRSeq_validations/validations_luciferase_plates_replicates_grid.csv", 
             header = T, 
             colClasses = "character", 
             na.strings = "", 
             fill= T)
rep <- melt(rep, id.vars = "rownames")
colnames(rep) <- c("row", "col", "tech_replicate")
scheme <- merge(scheme, rep)
scheme <- scheme[, .SD[, .(row, col, 
                           replicate= .GRP), 
                       .(tech_replicate,
                         date, 
                         plate)], 
                 Sample_ID]
scheme <- scheme[, !"tech_replicate"]

# ID/sample correspondance
pl <- as.data.table(read_xlsx("../../exp_data/vl_plasmids.xlsx"))[Experiment=="peSTARRSeq_validations"]
pl[, Sample_ID:= gsub(".* (.*)$", "\\1", labbook)]
scheme[pl, c("L", "R"):= tstrsplit(i.contains, "__SCR1__"), on= "Sample_ID"]

#------------------------------------------------------------#
# 2- Import luciferase data
#------------------------------------------------------------#
luc <- data.table(file= list.files("db/luciferase/peSTARRSeq_validations/", 
                                   pattern = "peSTARRvalid", 
                                   recursive = T, 
                                   full.names = T))
luc[, c("date", "plate", "lum") := tstrsplit(basename(file), "_|[.]", keep= c(1,3,4)), file]
luc <- luc[, melt(fread(file, header = T)[, 1:25], id.vars = "V1"), (luc)]
colnames(luc)[5:6] <- c("row", "col")
luc <- dcast(luc, date+plate+row+col~lum)
colnames(luc)[5:6] <- c("luc", "ren")

#------------------------------------------------------------#
# 3- Merge and process
#------------------------------------------------------------#
merged <- merge(scheme, luc, by= c("date", "plate", "row", "col"))
# Cutoffs rennilla and N tech replicates
dat <- merged[ren>7500]
dat[, check := .N>=3 & !is.na(Sample_ID), Sample_ID]
dat <- dat[(check), !"check"]
# Mean replicates
dat <- dat[, .(luc_norm= mean(luc/ren)), .(Sample_ID, L, R, bio_replicate= replicate)]

#------------------------------------------------------------#
# 4- Process similarly to STARR-Seq norm and merge with STARR-Seq
#------------------------------------------------------------#
# log2FC and sd
final <- dat[, .(log2FoldChange_luc= mean(log2(luc_norm), na.rm= T), 
                 sd_luc= sd(log2(luc_norm), na.rm= T)), .(L, R)]
# Center on negative controls
center <- mean(final[grepl("control", L) & grepl("control", R), log2FoldChange_luc])
final[, log2FoldChange_luc:= log2FoldChange_luc-center]
final[, log2FoldChange_luc_all:= .(.(dat[.BY, log2(luc_norm), on=c("L", "R")]-center)), .(L, R)]
# add features
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
final <- feat$add_feature(final, feat$lib)
# Add STARR
STARR <- fread("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.txt")
final <- merge(final, STARR, by= c("L", "R"), all.x = T)
# Add leftright/add act for barplot
final[, mean_luc_L:= mean(log2FoldChange_luc[group_R=="control"], na.rm= T), L]
final[, mean_luc_R:= mean(log2FoldChange_luc[group_L=="control"], na.rm= T), R]
final[, additive_luc:= log2(2^mean_luc_L+2^mean_luc_R)]

saveRDS(final, "Rdata/validations_luciferase_final_table.rds")
