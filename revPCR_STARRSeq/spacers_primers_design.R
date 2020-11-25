setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(Biostrings)

# Import genome sequences
seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)
seq <- seq[names(seq) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")]

# Make primers object
name <- c("revSPA1_300", "revSPA1_2k", "revSPA1_5k",
          "revSPA2_300", "revSPA2_2k", "revSPA2_5k",
          "revSPA3_300", "revSPA3_2k", "revSPA3_5k")
name <- paste0(rep(name, each= 2), c= c("_F", "_R"))
primers <- data.table(name= name, 
                      primer_seq= c("GCAACGATTTAATTAGGCGCTG", "TAGCTTGGCATGTTGATTCCAG",
                                    "CAGGCTTTCAACAAGTGCATTC", "GTCACAACTTTTCACACGCTG",
                                    "CCATCCATGCTCATAACCATCC", "CAACAGGAGCAACAACAACAAC",
                                    "TCAATTGATTTTAATGGCCCGC", "CGGAGATGGAAGATGGGAGTTC",
                                    "AGGCTATATCCCACGGTACATG", "CAAAGCCCGGGTATTATGTGAG",
                                    "CGGGATTGCGATTTCAAATTCC", "ATTGTGTCAGCTTAGTACGTGC",
                                    "AAGTCACTGGCTGCTGTAAATG", "TTTTCAGTCCGCTTATATGGCC",
                                    "CGGGTTTTCTTCCTTGCAATTG", "GTGTCGGATCTCTGTCAACATC",
                                    "CACAACTTTGATGGATTGCACG", "TTTATAGTTGTGGCTGCAGTGG"))
primers[, .id:= paste0("PCR", .GRP), tstrsplit(name, "_", keep=2)]
primers[, pattern:= ifelse(grepl("_F$", name), primer_seq, as.character(reverseComplement(DNAString(primer_seq)))), primers]
fwrite(primers, "Rdata/revPCR_STARRSeq/spacers_PCR_primers.txt", col.names = T, row.names = F, sep= "\t", quote= F)

# match NotI and primers patterns
.m <- rbind(data.table(name= "NotI", pattern= "GCGGCCGC", .id= "NotI"), primers, fill= T)
peaks <- .m[, {
  .c <- rbindlist(lapply(vmatchPattern(pattern, seq), as.data.table), idcol = T)
  colnames(.c)[1] <- "seqnames"
  .c
}, .m]

.tracks <- append(lapply(split(peaks, peaks$.id), GRanges),
                  list("RNA-seq"= "../available_data/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                       "ATAC-Seq"= "../available_data/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                       DSCP_200bp= "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                       DSCP_600bp= "../gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                       RPS12_200bp= "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                       RPS12_600bp= "../gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"))

spacers <- fread("Rdata/revPCR_STARRSeq/spacers_sequence.txt")

# previous script 34~40s
pdf("pdf/revPCR_STARRSeq/screenshot_spacers_PCRs.pdf", width = 25, height = 7)
par(mar= c(5,20,2,2))
spacers[, {
  .c <- GRanges(.SD)
  .c <- unlist(GRangesList(resize(.c, 1000, "center"), .c, resize(.c, 20000, "center"), resize(.c, 100000, "center")))
  my_screenshot(bw_GR_list = .tracks, bed = .c, gene.n = 5)
}, seq(nrow(spacers))]
dev.off()

