setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(GenomicRanges)

# Import previous lib (example)
link <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))

#-------------------------------#
# Mutant library (A)
#-------------------------------#
# Import mutant library
mut <- readRDS("db/library_design/twist015/mutant_sublib_twist015.rds")
mut[strand=="", c("seqnames", "strand"):= .(NA, NA)]
# Make and ID corresponding of group_mutation_A_sequenceNumber (sequence number shared between WT and mutated variants) 
mut[, ID:= paste0(group, "_", mut, "_A_", stringr::str_pad(as.numeric(as.factor(BAID)), 4, pad = "0")), group]
mut[, sublib:= "A"]
mut[, fw_linker:= "TTGACAGTGAGCGCGTCTCTCACCG"]
mut[, rev_linker:= "CCTAGGATCGACGCGGACAA"]
# Make enhancer start and end seq. unique (add AT if sequenced parts are similar)
mut[, uniq_start:= seq(.N), gsub("^(.{110}).*", "\\1", enh_sequence)]
mut[, uniq_end:= seq(.N), gsub("^.*(.{110})$", "\\1", enh_sequence)]
mut[, uniq_enh_sequence:= paste0(ifelse(uniq_start==2, "AT", ""), 
                                 enh_sequence,
                                 ifelse(uniq_end==2, "AT", "")), .(uniq_start, uniq_end)]
mut[, oligo_full_sequence:= paste0(fw_linker, uniq_enh_sequence, rev_linker)]
# Add random nt to reach 300
mut[nchar(oligo_full_sequence)==294, oligo_full_sequence:= paste0("GTA", oligo_full_sequence, "ATG")]
mut[nchar(oligo_full_sequence)==296, oligo_full_sequence:= paste0("GT", oligo_full_sequence, "TG")]
mut[nchar(oligo_full_sequence)==298, oligo_full_sequence:= paste0("G", oligo_full_sequence, "G")]

#-------------------------------#
# DHS library (B)
#-------------------------------#
DHS <- readRDS("db/library_design/twist015/DHS_sublib_twist015.rds")
DHS[, ID:= switch(type, 
                  "control"= "control",
                  "DHS"= "DHS",
                  "STARR"= "dev"), type]
DHS[, ID:= paste0(ID, "_B_", stringr::str_pad(1:1000, 4, pad = "0"))]
DHS[, sublib:= "B"]
DHS[, fw_linker:= "TTGACAGTGAGCGCGTCTCTCACCG"]
DHS[, rev_linker:= "GCTTTTGAAGCGTGCAGAATGAA"]
DHS[, oligo_full_sequence:= paste0("G", fw_linker, enh_sequence, rev_linker, "TA")]

# Assemble lib
lib <- rbind(mut[, .(sublib, ID, BAID, sublib, seqnames, start, end, strand, enh_sequence, oligo_full_sequence)],
             DHS[, .(sublib, ID, sublib, seqnames, start, end, strand, enh_sequence, oligo_full_sequence)],
             fill= T)
table(lib$strand, useNA= "ifany")
table(nchar(lib$oligo_full_sequence), useNA= "ifany")
length(unique(lib$ID))==nrow(lib)
lib[, .(length(unique(gsub({"^(.{140}).*"}, "\\1", oligo_full_sequence))), .N), sublib]
lib[, .(length(unique(gsub({"^.*(.{140})$"}, "\\1", oligo_full_sequence))), .N), sublib]
table(nchar(lib$oligo_full_sequence))

saveRDS(lib, 
        "Rdata/vl_library_twist015_112022.rds")
fwrite(lib[, .(ID, oligo_full_sequence)],
       "db/for_Bernado/vl_library_twist015_112022.txt", 
       col.names = T,
       sep= "\t",
       quote= F,
       na= NA)
