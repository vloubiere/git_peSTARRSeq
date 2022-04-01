setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(readxl)

#------------------------------------------------------------#
# 1- Import and format data
#------------------------------------------------------------#
dat <- data.table(Sample_ID_file= list.files("db/luciferase/peSTARRSeq_validations/", "scheme", recursive = T, full.names = T))
dat[, c("date", "plate") := tstrsplit(basename(Sample_ID_file), "_|[.]", keep= c(1,3)), Sample_ID_file]
dat[, tech_rep_file:= "db/luciferase/peSTARRSeq_validations/validations_luciferase_plates_replicates_grid.csv"]
dat[, luc_file:= list.files("db/luciferase/peSTARRSeq_validations/", 
                            pattern = paste0(date, ".*", plate, ".*lum1.csv$"),
                            recursive = T, 
                            full.names = T), .(date, plate)]
dat[, ren_file:= list.files("db/luciferase/peSTARRSeq_validations/", 
                            pattern = paste0(date, ".*", plate, ".*lum2.csv$"),
                            recursive = T, 
                            full.names = T), .(date, plate)]
cols <- grep("file$", names(dat), value= T)
dat[, gsub("_file$", "", cols):= lapply(.SD, function(x){
  .c <- fread(x, 
              header = T, 
              colClasses = "character",
              na.strings = as.character(NA), 
              fill= T)
  .c <- .c[, 1:25, with= F]
  setnames(.c, 1, "row")
  cols <- names(.c)[-1]
  if(!grepl("_scheme_", x))
    .c[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]
  .(melt(.c, id.vars = "row")[[3]])
}), (dat), .SDcols= cols]
dat <- dat[, lapply(.SD, unlist), .(plate, date), .SDcols= c("Sample_ID", "tech_rep", "luc", "ren")]

#------------------------------------------------------------#
# 2- Compute enrichment
#------------------------------------------------------------#
dat <- dat[ren>2500]
dat[, N_tech_rep:= .N, .(Sample_ID, plate, date)]
dat <- dat[N_tech_rep>2]
dat <- dat[, .(norm= mean(luc/ren)), .(Sample_ID, plate, date)]
dat[, N_bio_rep:= .N, Sample_ID]
dat <- dat[N_bio_rep>2, .(Sample_ID, norm)]

#------------------------------------------------------------#
# 3- Final table processed similarly to peSTARRSeq
#------------------------------------------------------------#
# ID/sample correspondance
samples <- as.data.table(read_xlsx("../../exp_data/vl_plasmids.xlsx"))[Experiment=="peSTARRSeq_validations"]
samples[, Sample_ID:= gsub(".* (.*)$", "\\1", labbook)]
dat[samples, c("L", "R"):= tstrsplit(i.contains, "__SCR1__"), on= "Sample_ID"]
# Log2 FoldChange
dat <- dat[, .(log2FoldChange= mean(log2(norm)),
               sd= sd(log2(norm))), .(L, R)]
# center on negative controls
dat[, log2FoldChange:= log2FoldChange-mean(dat[grepl("control", L) & grepl("control", R), log2FoldChange])]
# Compute individual act
dat[, mean_L:= mean(log2FoldChange[grepl("control", R)]), L]
dat[, mean_R:= mean(log2FoldChange[grepl("control", L)]), R]
dat[, additive:= log2(2^mean_L+2^mean_R)]
dat <- na.omit(dat)

saveRDS(dat, "Rdata/validations_luciferase_final_table.rds")
