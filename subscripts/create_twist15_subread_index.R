twist15 <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist015_112022.rds")
mutant <- twist15[sublib=="A"]
DHS <- twist15[sublib=="A"]

# Mutant
seqs <- mutant$oligo_full_sequence
names <- mutant$ID
seqinr::write.fasta(as.list(seqs), 
                    names, 
                    file.out = "db/fasta/twist15_mutant_lib.fasta",
                    as.string = T)
dir <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_mutant_lib/"
dir.create(dir, 
           showWarnings = F)
Rsubread::buildindex(basename = paste0(dir, "twist15"), 
                     reference = "db/fasta/twist15_mutant_lib.fasta")

# DHS
seqs <- DHS$oligo_full_sequence
names <- DHS$ID
seqinr::write.fasta(as.list(seqs), 
                    names, 
                    file.out = "db/fasta/twist15_DHS_lib.fasta",
                    as.string = T)
dir <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_DHS_lib/"
dir.create(dir, 
           showWarnings = F)
Rsubread::buildindex(basename = paste0(dir, "twist15"), 
                     reference = "db/fasta/twist15_DHS_lib.fasta")
