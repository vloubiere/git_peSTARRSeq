setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(readxl)

# Retrieve files ----
dat <- data.table(Sample_ID_file= list.files("db/luciferase/peSTARRSeq_validations/", "scheme", recursive = T, full.names = T))
dat[, c("date", "plate") := tstrsplit(basename(Sample_ID_file), "_|[.]", keep= c(1,3)), Sample_ID_file]
dat[, luc_file:= list.files("db/luciferase/peSTARRSeq_validations/", 
                            pattern = paste0(date, ".*", plate, ".*lum1.csv$"),
                            recursive = T, 
                            full.names = T), .(date, plate)]
dat[, ren_file:= list.files("db/luciferase/peSTARRSeq_validations/", 
                            pattern = paste0(date, ".*", plate, ".*lum2.csv$"),
                            recursive = T, 
                            full.names = T), .(date, plate)]
dat <- melt(dat, id.vars = c("date", "plate"))

# Import data ----
dat <- dat[, {
  # Import the different files
  .i <- .SD[,{
    .c <- fread(value, 
                header = T, 
                colClasses = "character",
                na.strings = as.character(NA), 
                sel= 1:25,
                fill= T)
    setnames(.c, 1, "row")
    melt(.c, id.vars = "row", variable.name= "column")
  }, variable]
  # dcast per row/column
  .i <- dcast(.i, row+column~variable, value.var = "value")
  setnames(.i, function(x) gsub("_file$", "", x))
  # Biological replicate
  .i$rep <- .GRP
  .i
}, .(date, plate)]
dat[, luc:= as.numeric(luc)]
dat[, ren:= as.numeric(ren)]

# Add enhancer IDs ----
samples <- as.data.table(read_xlsx("Rdata/vl_plasmids.xlsx"))[Experiment=="peSTARRSeq_validations"]
samples[, Sample_ID:= gsub(".* (.*)$", "\\1", labbook)]
dat[samples, c("L", "R"):= tstrsplit(i.contains, "__SCR1__"), on= "Sample_ID"]
# dat$Sample_ID <- NULL

# Remove points for which not enough replicates/signal ----
dat <- dat[ren>2500]
dat[, N_tech_rep:= .N, .(Sample_ID, rep)]
dat <- dat[N_tech_rep>2]
dat[, N_bio_rep:= length(unique(rep)), .(L, R)]
dat <- dat[N_bio_rep>2]
dat$N_bio_rep <- dat$N_tech_rep <- NULL

# Compute activity ----
dat <- dat[, .(log2FoldChange= mean(log2(luc/ren))), .(L, R, rep)]
dat <- dat[, .(log2FoldChange= mean(log2FoldChange),
               sd= sd(log2FoldChange)), .(L, R)]
# Center on negative controls
dat[, log2FoldChange:= log2FoldChange-mean(dat[grepl("^control", L) & grepl("^control", R), log2FoldChange])]

# Save ----
saveRDS(dat, "Rdata/validations_luciferase_final_table.rds")
