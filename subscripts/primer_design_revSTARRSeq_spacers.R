setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(Biostrings)

# Import genome sequences
seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)
seq <- seq[names(seq) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")]

#---------------------------------------#
# Make primers object
#---------------------------------------#
primers <- fread("db/library_design/alternative_spacers/selected_spacers_sequences.txt", 
                 select = c(1,2,3,8,9,12,10))
colnames(primers)[c(1,2,3,7)] <- paste0("region_", colnames(primers)[c(1,2,3,7)])
primers <- primers[rep(seq(nrow(primers)), each= 3)]
primers[, id:= c("revSPA1_300", "revSPA1_2k", "revSPA1_5k",
                 "revSPA2_300", "revSPA2_2k", "revSPA2_5k",
                 "revSPA3_300", "revSPA3_2k", "revSPA3_5k")]
primers <- primers[, .(name= paste0(id, c("_F", "_R"))), (primers)]
primers[, primer_seq:= c("GCAACGATTTAATTAGGCGCTG", "TAGCTTGGCATGTTGATTCCAG",
                         "CAGGCTTTCAACAAGTGCATTC", "GTCACAACTTTTCACACGCTG",
                         "CCATCCATGCTCATAACCATCC", "CAACAGGAGCAACAACAACAAC",
                         "TCAATTGATTTTAATGGCCCGC", "CGGAGATGGAAGATGGGAGTTC",
                         "AGGCTATATCCCACGGTACATG", "CAAAGCCCGGGTATTATGTGAG",
                         "CGGGATTGCGATTTCAAATTCC", "ATTGTGTCAGCTTAGTACGTGC",
                         "AAGTCACTGGCTGCTGTAAATG", "TTTTCAGTCCGCTTATATGGCC",
                         "CGGGTTTTCTTCCTTGCAATTG", "GTGTCGGATCTCTGTCAACATC",
                         "CACAACTTTGATGGATTGCACG", "TTTATAGTTGTGGCTGCAGTGG")]
primers[, primer_pattern:= ifelse(grepl("_F$", name), primer_seq, as.character(reverseComplement(DNAString(primer_seq)))), primers]

# Primers genomic location
primers <- rbindlist(lapply(seq(primers$primer_pattern), function(i) {
  res <- as.data.table(unlist(vmatchPattern(primers$primer_pattern[i], seq)))[, .(primer_seqnames= names, primer_start= start, primer_end= end)]
  res <- data.table(primers[i], res)
  return(res)
}))

# Amplicon sequence
primers[, amplicon_seq:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3, primer_seqnames, min(primer_start), max(primer_end))), .(primer_seqnames, id)]

fwrite(primers,
       "db/library_design/alternative_spacers/selected_spacers_PCR_primers.txt",
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F)
