setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir.create("db/sanger_sequencing/vllib015_016/", showWarnings = F)
require(data.table)
require(sangerseqR)
require(vlfunctions)
require(stringr)

# Result sumup
# 5 colonies sequences for vllib015 (DSCP) -> 3 OK, 2 unique enhancer
# 5 colonies sequences for vllib016 (RpS12) -> 2 OK, 2 unique enhancer, 1 bad sequencing

#---------------------------------------------------#
# IMPORT FILES and RENAME
# LEO did it for me -> copy files from his folder and rename them
#---------------------------------------------------#
dir.create("db/sanger_sequencing/vllib015_016_20210630/", 
           showWarnings = F)
files <- list.files("db/sanger_sequencing/vllib015_016_20210630/",
                    full.names = T)
file.remove(files)
files <- list.files("/groups/stark/serebreni/Vincent_colony_sanger/",
                    full.names = T)
file.copy(files, 
          paste0("db/sanger_sequencing/vllib015_016_20210630/", basename(files)))
files <- list.files("db/sanger_sequencing/vllib015_016_20210630/",
                    full.names = T)
for(i in 1:10)
{
  file.rename(files, 
              gsub(paste0("Colony", i, "_"), 
                   paste0("Colony", formatC(i, width = 2, flag = "00"), ifelse(i<6, "_vllib015_", "_vllib016_")), 
                   files))
  files <- list.files("db/sanger_sequencing/vllib015_016_20210630/",
                      full.names = T)
}

#---------------------------------------------------#
# IMPORT 
#---------------------------------------------------#
constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
lib <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))
dat <- data.table(file= list.files("db/sanger_sequencing/vllib015_016_20210630/", 
                                   full.names = T, 
                                   recursive = T))
dat[, ID:= gsub("SERE_(.*)_.*_.*$", "\\1", basename(file))]

#---------------------------------------------------#
# Identify enhancer sequences 
#---------------------------------------------------#
dat[, seq_L:= as.character(primarySeq(readsangerseq(grep("CAseq044", file, value = T)))), ID]
dat[, pat_L:= substr(seq_L, 
                     str_locate(seq_L, constructs["Flink_-2", sequence])[2]+1,
                     str_locate(seq_L, constructs["Flink_-2", sequence])[2]+21), seq_L]
dat[, seq_R:= as.character(primarySeq(readsangerseq(grep("CAseq003", file, value = T)))), ID]
dat[, pat_R:= substr(seq_R, 
                     str_locate(seq_R, vl_revComp(constructs["R1link+0", sequence]))[2]+1,
                     str_locate(seq_R, vl_revComp(constructs["R1link+0", sequence]))[2]+21), seq_R]
dat[is.na(pat_R), pat_R:= substr(seq_R, 
                                 str_locate(seq_R, vl_revComp(constructs["R2link+0", sequence]))[2]+1,
                                 str_locate(seq_R, vl_revComp(constructs["R2link+0", sequence]))[2]+21), seq_R]
dat[is.na(pat_R), pat_R:= substr(seq_R, 
                                 str_locate(seq_R, vl_revComp(constructs["R3link-3", sequence]))[2]+1,
                                 str_locate(seq_R, vl_revComp(constructs["R3link-3", sequence]))[2]+21), seq_R]
dat[!is.na(pat_L), enh_L:= lib[grep(paste0("^", pat_L), enh_seq), oligo_full_sequence], pat_L]
dat[is.na(pat_L), enh_L:= lib[1, oligo_full_sequence]]
dat[!is.na(pat_R), enh_R:= lib[grep(paste0(vl_revComp(pat_R), "$"), enh_seq), oligo_full_sequence], pat_R]
dat[is.na(enh_R), enh_R:= lib[1, oligo_full_sequence]]

#-------------------------------------------------------#
# Assemble full plasmid sequence
#-------------------------------------------------------#
DSCP_seq <- vl_digest(constructs["DSCP_STARRSeq", sequence], c("BshTI", "SalI"))
RPS_seq <- vl_digest(constructs["RPS12_STARRSeq", sequence], c("BshTI", "SalI"))
dat[, refseq:= paste0(ifelse(grepl("vllib015", ID), DSCP_seq[1], RPS_seq[1]), 
                      constructs["illumina_F", sequence],
                      enh_L, 
                      constructs["CGCov_F", sequence],
                      constructs["SCR1", sequence],
                      constructs["CGCov_R", sequence],
                      enh_R, 
                      constructs["illumina_R", sequence],
                      ifelse(grepl("vllib015", ID), DSCP_seq[3], RPS_seq[3]), 
                      collapse= ""), .(enh_L, enh_R, ID)]
dat[, refseq:= paste0(substr(refseq, 4001, nchar(refseq)),
                      substr(refseq, 1, 4000)), refseq]
dat[, rev:= ifelse(grepl("CAseq001|CAseq044", file), F, T)]

#-------------------------------------------------------#
# Final check and plot
#-------------------------------------------------------#
pdf("pdf/sanger_sequencing/sanger_vllib015_016_20210630.pdf")
dat[, vl_sanger_align(refseq, 
                      abfiles = file, 
                      revcomp = rev, 
                      feat_sequences = constructs[c("DSCP", "Rps12", "SCR1"), sequence]), refseq]
dev.off()
