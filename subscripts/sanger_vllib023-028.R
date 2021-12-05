setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir <- "db/sanger_sequencing/vllib023_028/"
dir.create(dir, showWarnings = F)
require(vlfunctions)
require(sangerseqR)
require(stringr)

#---------------------------------------------------#
# IMPORT FILES AND COPY THEM LOCALLY
#---------------------------------------------------#
# source("../../exp_data/update_files.R")
dat <- as.data.table(read_xlsx(path = "/groups/stark/vloubiere/exp_data/vl_sanger_sequencing.xlsx"))
dat <- dat[grepl("vllib", name) & Date==211203 & project=="pe_STARRSeq"]
dat[, molBio:= list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-12-03_PCW-IMP-253_0701/",
                          paste0("_", seqID, "_"),
                          full.names = T), seqID]
dat[, file:= list.files(dir,
                        paste0("_", seqID, "_"),
                        full.names = T), seqID]
dat[is.na(file) & file.exists(molBio), file:= paste0(dir, basename(molBio))]
dat[!file.exists(file), file.copy(molBio, file), .(molBio, file)]

#---------------------------------------------------#
# IMPORT sequences
#---------------------------------------------------#
constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
pl <- as.data.table(read_excel("../../exp_data/vl_plasmids.xlsx"), 
                    key= "ID")
lib <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))

#---------------------------------------------------#
# Identify enhancer sequences 
#---------------------------------------------------#
dat[, seq_L:= as.character(primarySeq(readsangerseq(grep("CASeq044", file, value = T)))), Sample_ID]
dat[, pat_L:= substr(seq_L, 
                     str_locate(seq_L, constructs["Flink_-2", sequence])[2]+1,
                     str_locate(seq_L, constructs["Flink_-2", sequence])[2]+21), seq_L]
dat[, seq_R:= as.character(primarySeq(readsangerseq(grep("CASeq003", file, value = T)))), Sample_ID]
dat[, pat_R:= substr(seq_R, 
                     str_locate(seq_R, constructs["R1link+0", sequence])[2]+1,
                     str_locate(seq_R, constructs["R1link+0", sequence])[2]+21), seq_R]
dat[is.na(pat_R), pat_R:= substr(seq_R, 
                                 str_locate(seq_R, constructs["R2link+0", sequence])[2]+1,
                                 str_locate(seq_R, constructs["R2link+0", sequence])[2]+21), seq_R]
dat[is.na(pat_R), pat_R:= substr(seq_R, 
                                 str_locate(seq_R, constructs["R3link-3", sequence])[2]+1,
                                 str_locate(seq_R, constructs["R3link-3", sequence])[2]+21), seq_R]
dat[!is.na(pat_L), enh_L:= lib[grep(paste0("^", pat_L), enh_seq), oligo_full_sequence], pat_L]
dat[is.na(pat_L), enh_L:= lib[1, oligo_full_sequence]]
dat[!is.na(pat_R), enh_R:= lib[grep(paste0(vl_revComp(pat_R), "$"), enh_seq), oligo_full_sequence], pat_R]
dat[is.na(enh_R), enh_R:= lib[1, oligo_full_sequence]]

#-------------------------------------------------------#
# Assemble full plasmid sequence
#-------------------------------------------------------#
dat[, c("left" , "right"):= {
  .c <- vl_digest(pl[backbone, Sequence], c("BshTI", "SalI"))
  print(length(.c))
  .(.c[1], .c[3])
}, backbone]
dat[, refseq:= paste0(left,
                      constructs["illumina_F", sequence],
                      enh_L, 
                      constructs["CGCov_F", sequence],
                      ifelse(grepl("vllib023|vllib024", name), 
                             paste0("ATACGCGCGCG", constructs["intron4_shortened", sequence], "CGCGCGCGTAT"), 
                             constructs["SCR1", sequence]),
                      constructs["CGCov_R", sequence],
                      enh_R, 
                      vl_revComp(constructs["illumina_R", sequence]),
                      right, 
                      collapse= ""), .(enh_L, enh_R, left, right, Sample_ID, name)]
dat[, nchar(refseq), refseq]
dat[, refseq:= paste0(substr(refseq, nchar(refseq)-1000, nchar(refseq)), 
                      substr(refseq, 1, nchar(refseq)-1001)), refseq]

#-------------------------------------------# Final check and plot
#-------------------------------------------------------#
pdf("pdf/sanger_sequencing/sanger_vllib023-028_031221.pdf")
dat[, {
  vl_sanger_align(refseq, 
                  abfiles = file, 
                  revcomp = as.logical(revComp), 
                  feat_sequences = constructs[c("DSCP", "RPS12", "hk_high_2", "hk_low_1", "dev_high_1", "dev_low_FL8", "SCR1", "intron4_shortened"), sequence], 
                  feat_cols = c("blue", "red", "green", "gold", "cornflowerblue", "firebrick3", "black", "grey"))
  mtext(name, line = 3)
}, .(refseq, name)]
dev.off()
