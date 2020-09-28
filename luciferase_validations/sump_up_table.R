setwd("/groups/stark/vloubiere/projects/0006_vllib002_SCR1_validations_luc/")
source("../../scripts/R_functions/my_SANGER_aligner.R")
require(data.table)
require(colorRamps)
require(Biostrings)
require(sangerseqR)
require(readxl)
require(gridExtra)

dat <- as.data.table(read_xlsx("metadata/clean_stocks.xlsx"))
dat <- dat[seq_ID_1 != "NA", Sample_ID:enh_R]
dat <- dat[order(Sample_ID)]
dat[grepl("^ctl", ID_L) & grepl("^ctl", ID_R), type:= "control"]
dat[grepl("^ctl", ID_L) & !grepl("^ctl", ID_R), type:= "enh_R"]
dat[!grepl("^ctl", ID_L) & grepl("^ctl", ID_R), type:= "enh_L"]
dat[!grepl("^ctl", ID_L) & !grepl("^ctl", ID_R), type:= "enh_pair"]
dat <- dat[, .(Sample_ID, type, enh_L, enh_R, ID_L, ID_R)]
dat <- dat[, .(Sample_ID= .(Sample_ID)), .(enh_L, enh_R, type)]

dat[type=="control", use:= T]
dat[type=="enh_pair" & enh_L %in% dat[type=="enh_L", enh_L] & enh_R %in% dat[type=="enh_R", enh_R], use:= T]
dat[type=="enh_L" & enh_L %in% dat[type=="enh_pair" & use=="TRUE", enh_L], use:= T]
dat[type=="enh_R" & enh_R %in% dat[type=="enh_pair" & use=="TRUE", enh_R], use:= T]
dat[is.na(use), use:= F]

saveRDS(dat, "Rdata/sum_up.rds")
fwrite(dat, "Rdata/sum_up.txt", col.names = T, row.names = F, sep= "\t", quote= F)

check <- dat[use=="TRUE" & lengths(Sample_ID)>1]
