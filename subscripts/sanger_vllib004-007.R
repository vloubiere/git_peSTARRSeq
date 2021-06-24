setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(seqinr)
require(gridExtra)
require(Biostrings)
require(sangerseqR)

constructs <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", 
                    key= "name")
enh_sequences <- read.fasta("db/fasta/vllib001-014.fasta", as.string = T)

#---------------------------------------#
# Retrieve all sequences
#---------------------------------------#
dat <- data.table(file= list.files("db/sanger_sequencing/revPCR_pilot", full.names = T))
dat[, ID:= tstrsplit(basename(file), "_", keep= 2)]
dat <- dat[!ID %like% "1.1.2"] # shitty point
# Variable sequences
L_ID <- c("0873$","0781$","0991$","0035$","0089$","0950$","0521$","0205$", rep("SCR2", 8))
dat[, L_name:= sapply(L_ID, function(x) grep(x, names(enh_sequences), value = T))[.GRP], ID]
dat[, L:= as.character(enh_sequences[L_name]), L_name]
R_ID <- c("0830$","0452$","0452$","0220$","0203$","0169$","0736$","0112$","0953$","0239$","0616$","0274$","0401$","0578$","0117$","0650$")
dat[, R_name:= sapply(R_ID, function(x) grep(x, names(enh_sequences), value = T))[.GRP], ID]
dat[, R:= as.character(enh_sequences[R_name]), R_name]
dat[, spacer:= c(rep(constructs["revSPA1_300", sequence], 4), rep(constructs["revSPA1_2k", sequence], 12))[.GRP], ID]
# Cat qith constant sequences
DSCP_dig <- vl_digest(constructs["DSCP_STARRSeq", sequence], c("BshTI", "SalI"))
dat[, refseq:= paste0(DSCP_dig[1], 
                      "CGCGCGCG", 
                      # constructs["Flink_-2", sequence],
                      L,
                      # ifelse(grepl("_B_", L_name), constructs["R2link+0", sequence], constructs["R1link+0", sequence]),
                      constructs["rev_illumina_F", sequence],
                      constructs["CGCov_F", sequence],
                      spacer,
                      constructs["CGCov_R", sequence],
                      constructs["rev_illumina_R", sequence],
                      # constructs["Flink_-2", sequence],
                      R,
                      # ifelse(grepl("_B_", L_name), constructs["R2link+0", sequence], constructs["R3link-3", sequence]),
                      DSCP_dig[3]), .(L_name, R_name)]


pdf("pdf/sanger_sequencing/revPCR_STARRSeq_pilot_constructs.pdf", 20, 10)
par(mfrow=c(2,2))
dat[, {
  vl_sanger_align(refseq = refseq, 
                  abfiles = file,
                  revcomp = ifelse(grepl("CASeq044", file), F, T),
                  feat_sequences = c(DSCP_dig[1], L, spacer, R, DSCP_dig[3]),
                  feat_names = c("upstream", "Left", "Spacer","Right", "downstream"))
}, .(refseq, L, spacer, R)]
dev.off()


# # MaubI
# DSCP_dig <- vl_digest(constructs["DSCP_STARRSeq", sequence], c("BshTI", "SalI"))
# upstream <- paste0(DSCP_dig[1], "CGCGCGCG", constructs["linker_F-2_TWIST", sequence])
# left_linker <- constructs[c("rev_illumina_F", "CGCov_F"), paste0(sequence, collapse = "")]
# right_linker <- constructs[c("CGCov_R", "rev_illumina_R", "linker_F-2_TWIST"), paste0(sequence, collapse = "")]
# downstream <- DSCP_dig[2]
# 
# # Variable regions
# lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
# spa <- fread("Rdata/selected_spacers_PCR_primers.txt")
# spa <- unique(spa[, .(SPA= id, SPA_seq= amplicon_seq)])
# spa[SPA %in% c("revSPA1_300", "revSPA1_2k", "revSPA3_2k")]
# 
# # Left enhancer
# dat <- rbind(lib[linker_ID=="B", .(L= ID_vl, L_seq= paste0(enh_sequence, constructs["R2link+0", sequence]))], 
#              data.table(L= "SCR2", L_seq= constructs[c("SCR2", "R2link+0"), paste0(sequence, collapse="")]))
# # SPACER
# dat <- dat[, (spa), (dat)]
# # Right enhancer
# dat <- dat[, {
#   if(L!="SCR2"){
#     lib[linker_ID=="B", .(R= ID_vl, R_seq= paste0(enh_sequence, constructs["R2link+0", sequence]))]
#   }else if(L=="SCR2"){
#     lib[linker_ID=="C", .(R= ID_vl, R_seq= paste0(enh_sequence, constructs["R3link-3", sequence]))]
#   }
# }, (dat)]
# dat[, construct_sequence:= paste0(upstream, "CGCGCGCG", constructs["Flink_-2", sequence], L_seq, left_linker, SPA_seq, right_linker, R_seq, downstream)]
# dat[, full_sequence:= paste0(upstream, L_seq, left_linker, SPA_seq, right_linker, R_seq, downstream)]
# 
# # DONT NEED TO BE RE_RUN, see identified sequences hardcoded below
# # #---------------------------------------#
# # # Match sanger results
# # #---------------------------------------#             NOT FULLY EFFECTIVE!!
# # # Score
# # if(!file.exists("Rdata/revPCR_STARRSeq/sanger_alignment_vllib004_5_6_7.rds")){
# #   files <- data.table(ab_F= list.files("db/sanger_sequencing/revPCR_lib004_5_6_7/", "vl251|vl253|vl263|vl265|CASeq044.*.ab1$", full.names = T),
# #                       ab_R= list.files("db/sanger_sequencing/revPCR_lib004_5_6_7/", "vl254|vl266|GSP1.*.ab1$", full.names = T))
# #   match <- dat[, (files), (dat)]
# #   # Filter only meaningful combinations
# #   match <- match[grepl("_1.|_2.", ab_F)]
# #   match[, F_score:= {
# #     print(.GRP)
# #     pairwiseAlignment(as.character(primarySeq(readsangerseq(ab_F))), L_seq, scoreOnly= T)
# #   }, .(L_seq, ab_F)]
# #   match[, R_score:= {
# #     print(.GRP)
# #     pairwiseAlignment(as.character(reverseComplement(primarySeq(readsangerseq(ab_R)))), R_seq, scoreOnly= T)
# #   }, .(R_seq, ab_R)]
# #   saveRDS(match, "Rdata/revPCR_STARRSeq/sanger_alignment_vllib004_5_6_7.rds")
# # }
# # match <- readRDS("Rdata/revPCR_STARRSeq/sanger_alignment_vllib004_5_6_7.rds")
# 
# #---------------------------------------#
# # Screenshots
# #---------------------------------------#
# Left <- sapply(c("0873$","0781$","0991$","0035$","0089$","0950$","0521$","0205$", rep("SCR2$", 8)), function(x) unique(grep(x, dat$L, value = T)))
# Spacer <- c(rep("revSPA1_300", 4), rep("revSPA1_2k", 12))
# Right <- sapply(c("0830$","0452$","0452$","0220$","0203$","0169$","0736$","0112$","0953$","0239$","0616$","0274$","0401$","0578$","0117$","0650$"), function(x) unique(grep(x, dat$R, value = T)))
# setkeyv(dat, c("L", "SPA", "R"))
# res <- dat[list(Left, Spacer, Right)]
# res[, ab_F:= list.files("db/sanger_sequencing/revPCR_pilot", "CASeq044.*.ab1$", full.names = T)]
# res[, ab_R:= list.files("db/sanger_sequencing/revPCR_pilot", "GSP1.*.ab1$", full.names = T)[-2]]
# 
# pdf("pdf/sanger_sequencing/revPCR_STARRSeq_pilot_constructs.pdf", 20, 10)
# par(mfrow=c(2,2))
# res[, {
#   vl_sanger_align(refseq = paste0(upstream, L_seq, left_linker, SPA_seq, right_linker, R_seq, downstream), 
#                   abfiles = c(ab_F, ab_R),
#                   revcomp = c(F, T),
#                   feat_sequences = c(upstream, L_seq, left_linker, SPA_seq, right_linker, R_seq, downstream),
#                   feat_names = c("upstream", "Left", "L linker", "Spacer", "R linker", "Right", "downstream"))
# }, (res)]
# dev.off()
