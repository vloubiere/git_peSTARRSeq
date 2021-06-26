require(data.table)

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
pl <- fread("../../exp_data/vl_plasmids.txt")[Experiment=="peSTARRSeq_validations"]
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
# Normalize for negative controls
dat[, luc_norm:= luc_norm/mean(dat[grepl("^control", L) & grepl("^control", R), luc_norm])]
# Compute additive scores
dat[, luc_mean_L:= mean(.SD[grepl("^control", R), luc_norm], na.rm= T), .(bio_replicate, L)]
dat[, luc_mean_R:= mean(.SD[grepl("^control", L), luc_norm], na.rm= T), .(bio_replicate, R)]
dat[, luc_add:= luc_mean_L+luc_mean_R]

saveRDS(dat, "Rdata/validations_luciferase_final_table.rds")