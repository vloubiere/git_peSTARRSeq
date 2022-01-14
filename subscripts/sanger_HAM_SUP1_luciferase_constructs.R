setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(plater)

constructs <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", key= "name")

pl <- as.data.table(readxl::read_xlsx("/groups/stark/vloubiere/exp_data/vl_plasmids.xlsx"))
pl <- pl[Experiment=="ham_pilot_luc"]

#----------------------------#
# Refseq
#----------------------------#
p1 <- "db/sanger_sequencing/ham_pilot/100220/"
pl[, file_1:= .(list(list.files(p1, x, full.names= T))), .(x= paste0("pluc", gsub(".* (.*$)", "\\1_", labbook)))]
p2 <- "db/sanger_sequencing/ham_pilot/240220_recloning/"
pl[, file_2:= .(list(list.files(p2, x, full.names = T))), .(x= paste0("_p", formatC(as.numeric(gsub(".* (.*$)", "\\1", labbook)), width = 2, flag = "0")))]
p3 <- "db/sanger_sequencing/ham_pilot/240220_reseq_last_sample/"
pl[, file_3:= .(list(list.files(p3, x, full.names= T))), .(x= paste0("pluc", gsub(".* (.*$)", "\\1_", labbook)))]
pl <- melt(pl, measure.vars = patterns("file"))
pl <- pl[lengths(value)>0, .(value= unlist(value)), setdiff(colnames(pl), "value")]
plasmids <- as.data.table(readxl::read_xlsx("/groups/stark/vloubiere/exp_data/vl_plasmids.xlsx"))
upstream <- vl_digest(plasmids[grepl("pLuc002", ID), Sequence], "KpnI")[1]
downstream <- vl_digest(plasmids[grepl("pLuc002", ID), Sequence], "KpnI")[2]
pl[grepl("-", contains), seq:= paste0(substr(upstream, 
                                             start = nchar(upstream)-1000,
                                             stop = nchar(upstream)),
                                      "TCGATACCGACACCATTGAG", # RTTAseqRmod zuber
                                      constructs["Flink_+0", sequence],
                                      constructs[gsub("(^.*)-(.*)-(.*$)", "\\1", contains), sequence],
                                      vl_revComp(constructs["R1link+0", sequence]),
                                      constructs["CGCov_F", sequence],
                                      constructs[gsub("(^.*)-(.*)-(.*$)", "\\2", contains), sequence],
                                      vl_revComp(constructs["CGCov_R", sequence]),
                                      constructs["Flink_+0", sequence],
                                      constructs[gsub("(^.*)-(.*)-(.*$)", "\\3", contains), sequence],
                                      vl_revComp(constructs["R1link+0", sequence]),
                                      "CGTAGATGTACTGCCAAGTAGTAC", # CMVseqFmod zuber
                                      substr(downstream, 1, 1000)), contains]
pl[!grepl("-", contains), seq:= paste0(substr(upstream, 
                                              start = nchar(upstream)-1000,
                                              stop = nchar(upstream)),
                                       "TCGATACCGACACCATTGAG", # RTTAseqRmod zuber
                                       constructs["Flink_+0", sequence],
                                       constructs[gsub("(^.*)-(.*)-(.*$)", "\\1", contains), sequence],
                                       vl_revComp(constructs["R1link+0", sequence]),
                                       "CGTAGATGTACTGCCAAGTAGTAC", # CMVseqFmod zuber
                                       substr(downstream, 1, 1000)), contains]
pl[, rev:= ifelse(grepl("CASeq001", value), F, T)]

#----------------------------#
# Plot
#----------------------------#
pdf("/groups/stark/vloubiere/projects/pe_STARRSeq/pdf/luciferase/sanger_sequencing.pdf", height = 5)
par(mar= c(2,20,5,2))
pl[, 
   {
     vl_sanger_align(refseq = seq, 
                     abfiles = value, 
                     revcomp = rev, 
                     feat_sequences = constructs[c("HAM1", "SCR1", "SUP1"), sequence], 
                     feat_names = c("HAM1", "SCR1", "SUP1"),
                     feat_cols = c("green", "black", "red"))
     mtext(contains)
   }, .(seq, dirname(value))]
dev.off()
