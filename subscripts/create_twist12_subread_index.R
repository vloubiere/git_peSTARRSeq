seqs <- twist12$oligo_full_sequence
names <- twist12$ID
write.fasta(as.list(seqs), 
            names, 
            file.out = "db/fasta/twist12_lib.fasta",
            as.string = T)
dir <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/"
dir.create(dir, 
           showWarnings = F)
Rsubread::buildindex(basename = paste0(dir, "twist12"), 
                     reference = "db/fasta/twist12_lib.fasta")