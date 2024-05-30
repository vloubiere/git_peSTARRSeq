setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import selected enhaners ----
dat <- fread("db/PCR_primers/PCR_primers_homotypic_pairs_luc.txt")[variable=="left_primer"]
dat[, number:= .I]

# sequences ----
lib <- readRDS("Rdata/vl_library_twist008_112019.rds")
lib <- as.data.table(lib)
dat[lib, seq:= oligo_full_sequence, on= "enh==ID_vl"]
