setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(readxl)

#------------------------------------------------------------#
# 1- Import and format data
#------------------------------------------------------------#
# Plate schemes
scheme <- data.table(file= list.files("db/luciferase/peSTARRSeq_validations/", "scheme", recursive = T, full.names = T))
scheme[, c("date", "plate") := tstrsplit(basename(file), "_|[.]", keep= c(1,3)), file]
scheme <- scheme[, melt(fread(file, header = T, colClasses = "character", na.strings = "", fill= T), 
                        id.vars = "rownames", value.name = "Sample_ID"), (scheme)]
colnames(scheme)[4:5] <- c("row", "col")
# Add technical replicates
rep <- fread("Rdata/luciferase_validations/plates_replicates_grid.csv", header = T, colClasses = "character", na.strings = "", fill= T)
rep <- melt(rep, id.vars = "rownames")
colnames(rep) <- c("row", "col", "tech_replicate")
scheme <- merge(scheme, rep)
scheme <- scheme[, .SD[, .(row, col, replicate= .GRP), .(tech_replicate, date, plate)], Sample_ID]
scheme <- scheme[, !"tech_replicate"]
# ID/sample correspondance
ID <- as.data.table(read_excel("Rdata/luciferase_validations/clean_stocks.xlsx"))
scheme[ID, c("enh_L", "enh_R"):= .(i.enh_L, i.enh_R), on= "Sample_ID"]
# luciferase data
luc <- data.table(file= list.files("db/luciferase/peSTARRSeq_validations/", "peSTARRvalid", recursive = T, full.names = T))
luc[, c("date", "plate", "lum") := tstrsplit(basename(file), "_|[.]", keep= c(1,3,4)), file]
luc <- luc[, melt(fread(file, header = T)[, 1:25], id.vars = "V1"), (luc)]
colnames(luc)[5:6] <- c("row", "col")
luc <- dcast(luc, date+plate+row+col~lum)
colnames(luc)[5:6] <- c("luc", "ren")

#------------------------------------------------------------#
# 2- Process luciferase data
#------------------------------------------------------------#
merged <- merge(scheme, luc, by= c("date", "plate", "row", "col"))
# Cutoffs rennilla and N tech replicates
dat_all <- merged[ren>7500]
dat_all[, check := .N>=3 & !is.na(Sample_ID), Sample_ID]
dat_all <- dat_all[(check), !"check"]
# Mean replicates
dat_all <- dat_all[, .(luc_norm= mean(luc/ren)), .(Sample_ID, enh_L, enh_R, replicate)]
# Normalize for negative controls
dat_all[, luc_norm:= luc_norm/mean(dat_all[grepl("^control", enh_L) & grepl("^control", enh_R), luc_norm])]
# Compute additive scores
dat_all[, luc_mean_L:= mean(.SD[grepl("^control", enh_R), luc_norm], na.rm= T), .(replicate, enh_L)]
dat_all[, luc_mean_R:= mean(.SD[grepl("^control", enh_L), luc_norm], na.rm= T), .(replicate, enh_R)]
dat_all[, luc_add:= luc_mean_L+luc_mean_R]
dat_all[grepl("^control", enh_L)|grepl("^control", enh_R), luc_add:= NA]

saveRDS(dat_all, "Rdata/luciferase_validations/C_luc_validations_final_table.rds")







