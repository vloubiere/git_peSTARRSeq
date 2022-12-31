twist15 <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist015_112022.rds")
seqs <- twist15$oligo_full_sequence
names <- twist15$ID
seqinr::write.fasta(as.list(seqs), 
                    names, 
                    file.out = "db/fasta/twist15_lib.fasta",
                    as.string = T)
dir <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_lib/"
dir.create(dir, 
           showWarnings = F)
Rsubread::buildindex(basename = paste0(dir, "twist15"), 
                     reference = "db/fasta/twist12_lib.fasta")
