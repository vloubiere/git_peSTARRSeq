setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import feat
feat <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
# Add chromosome sizes and compute sets
chr_sizes <- as.data.table(GenomicFeatures::getChromInfoFromUCSC("dm3"))
feat[chr_sizes, length:= i.length, on= "seqnames==chrom"]
feat[, set:= fcase(seqnames=="chr2R" & between(start, 1, length/2), "validation",
                   seqnames=="chr2R" & between(start, length/2, length), "test", 
                   default= "train"), .(seqnames, length)]

# Import vllib002 and merge with feat
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
lib[feat, c("set_L", "seq_L"):= .(i.set, i.enh_sequence), on= "L==ID_vl"]
lib[feat, c("set_R", "seq_R"):= .(i.set, i.enh_sequence), on= "R==ID_vl"]
# Keep only the lines for which both enhancers are in the same set
final <- lib[set_L==set_R, .(L, R, log2FoldChange, set= set_L, seq_L= seq_L, seq_R= seq_R)]

# Save
final[, seqinr::write.fasta(sequences = as.list(seq_L), 
                            names = L, 
                            file.out = paste0("db/for_Bernado/deepSTARR/", set, "_left.fasta"), 
                            as.string = T), set]
final[, seqinr::write.fasta(sequences = as.list(seq_R), 
                            names = R, 
                            file.out = paste0("db/for_Bernado/deepSTARR/", set, "_right.fasta"), 
                            as.string = T), set]
final[, fwrite(.(log2FoldChange),
               paste0("db/for_Bernado/deepSTARR/", set, "_activity.txt"),
               quote= F,
               col.names = F), set]
