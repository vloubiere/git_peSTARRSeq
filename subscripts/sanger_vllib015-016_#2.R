setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir.create("db/sanger_sequencing/vllib015_016/", showWarnings = F)
require(data.table)
require(sangerseqR)
require(vlfunctions)
require(stringr)
require(readxl)

# Result sumup
# vllib015 --> 4 cols are good, one has bad sequencing but also seems ok (benchling)
# vllib016 --> 2 cols are good, 2 have bad sequencing but also seem ok (benchling), 009 cannot say (sequencing is bad)

#---------------------------------------------------#
# IMPORT FILES 
#---------------------------------------------------#
dat <- as.data.table(read_xlsx(path = "/groups/stark/vloubiere/exp_data/vl_sanger_sequencing.xlsx"))
dat <- dat[grepl("vllib", name) & Date==210714 & project=="pe_STARRSeq"]
# files <- list.files(c("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-07-14_PCW-IMP-253_0231/",
#                       "/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-07-14_PCW-IMP-253_0232/"),
#                     full.names = T)
# files <- files[grepl(paste0(dat$seqID, collapse= "|"), files)]
# dir.create("db/sanger_sequencing/vllib015_016_new_20210714/",
#            showWarnings = F)
# file.copy(files,
#           paste0("db/sanger_sequencing/vllib015_016_new_20210714/", basename(files)))

#---------------------------------------------------#
# IMPORT sequences
#---------------------------------------------------#
# source("/groups/stark/vloubiere/exp_data/update_files.R")
constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
pl <- fread("../../exp_data/vl_plasmids.txt", 
            key= "ID")
lib <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))

#---------------------------------------------------#
# Identify enhancer sequences 
#---------------------------------------------------#
dat[, file:= list.files("db/sanger_sequencing/vllib015_016_new_20210714/", 
                        seqID, 
                        full.names = T), seqID]
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
DSCP_seq <- vl_digest(pl["STARRSeq_pGL3_IF_DSCPII_mhcI_GFP_AgeI-SalI", Sequence], c("BshTI", "SalI"))
RPS_seq <- vl_digest(pl["pGL3_IF_RpS12_mhcI_GFP_AgeI-SalI", Sequence], c("BshTI", "SalI"))
dat[, refseq:= paste0(ifelse(grepl("vllib15", name), DSCP_seq[1], RPS_seq[1]), 
                      constructs["illumina_F", sequence],
                      enh_L, 
                      constructs["CGCov_F", sequence],
                      constructs["SCR1", sequence],
                      constructs["CGCov_R", sequence],
                      enh_R, 
                      vl_revComp(constructs["illumina_R", sequence]),
                      ifelse(grepl("vllib15", name), DSCP_seq[3], RPS_seq[3]), 
                      collapse= ""), .(enh_L, enh_R, Sample_ID, name)]
dat[, refseq:= paste0(substr(refseq, 4001, nchar(refseq)),
                      substr(refseq, 1, 4000)), refseq]

#-------------------------------------------------------#
# Final check and plot
#-------------------------------------------------------#
pdf("pdf/sanger_sequencing/sanger_vllib015_016_new_20210714.pdf")
# par(mar= c(5,10,5,5))
dat[, vl_sanger_align(refseq, 
                      abfiles = file, 
                      revcomp = as.logical(revComp), 
                      feat_sequences = constructs[c("DSCP", "RPS12", "SCR1"), sequence]), refseq]
dev.off()
