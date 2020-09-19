setwd("/groups/stark/vloubiere/projects/0006_vllib002_SCR1_validations_luc/")
source("/groups/stark/vloubiere/scripts/R_functions/my_plots.R")
require(data.table)
require(Biostrings)

dat <- fread("Rdata/A1_selected_extra_sup_20200625.txt")
lib <- as.data.table(readRDS("../pe_STARRSeq/Rdata/A_vl_library_112019.rds"))

dat[lib, sequence:= i.enh_sequence, on= "enh_ID==ID_vl"]
dat[position=="left", primer:= substr(sequence, 1, 20)]
dat[position=="left", pluc_overhang:= "ATTTCTCTATCGATAGGTAC"]
dat[position=="right", primer:= as.character(reverseComplement(DNAStringSet(substr(sequence, nchar(sequence)-19, nchar(sequence)))))]
dat[position=="right", pluc_overhang:= "GACGCGTAAGAGCTCGGTAC"]

dat[, primer_overhang := paste0(pluc_overhang, primer)]
dat[position=="left", primer_name := paste0(valid_ID, "_", "F")]
dat[position=="right", primer_name := paste0(valid_ID, "_", "R")]

fwrite(dat, paste0("Rdata/B1_primers_to_order_", gsub("-", "", Sys.Date()), ".txt"))