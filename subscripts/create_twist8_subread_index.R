twist8 <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist008_112019.rds")
seqs <- c(twist8$oligo_full_sequence,
          paste0(constructs["Flink_lib", sequence], 
                 constructs[c("SCR1", "SCR2", "HAM1", "SUP1"), sequence], 
                 constructs["R1link_lib", sequence]))
names <- c(twist8$ID_vl,
           "ts_SCR1_01001", 
           "ts_SCR2_01002",
           "ts_HAM1_01003",
           "ts_SUP1_01004")
seqinr::write.fasta(as.list(seqs), 
                    names, 
                    file.out = "db/fasta/twist8_lib.fasta",
                    as.string = T)
dir <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/"
dir.create(dir, 
           showWarnings = F)
Rsubread::buildindex(basename = paste0(dir, "twist8"), 
                     reference = "db/fasta/twist8_lib.fasta")