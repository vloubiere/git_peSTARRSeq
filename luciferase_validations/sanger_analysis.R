setwd("/groups/stark/vloubiere/projects/0006_vllib002_SCR1_validations_luc/")
source("../../scripts/R_functions/my_SANGER_aligner.R")
require(data.table)
require(colorRamps)
require(Biostrings)
require(sangerseqR)
require(readxl)

# Import candidates
cand <- data.table(file= c("Rdata/A_selected_candidates_20200521.txt", "Rdata/A1_selected_extra_sup_20200625.txt"))
cand <- setDT(cand[, fread(file), file], key = "position")

# Assemble refseqs
lib <- as.data.table(readRDS("../pe_STARRSeq/Rdata/vl_library_112019.rds"))
constructs <- fread("../../data/constructs_sequences/constructs_sequences.txt", key= "name")

seq <- CJ(ID_L= cand["left", valid_ID], ID_R= cand["right", valid_ID])
seq[cand["left"], enh_L:= i.enh_ID, on= "ID_L==valid_ID"]
seq[cand["right"], enh_R:= i.enh_ID, on= "ID_R==valid_ID"]
seq[lib, seq_L:= i.enh_sequence, on= "enh_L==ID_vl"]
seq[lib, seq_R:= i.enh_sequence, on= "enh_R==ID_vl"]
seq[grepl("_A_", enh_L), refseq:= paste0(constructs["CASeq001_to_F_gibson_ov", sequence], seq_L, 
                                         constructs["rev_comp_TWIST_linker_R1", sequence],
                                         constructs["SCR1_ov", sequence], 
                                         constructs["TWIST_linker_F", sequence], seq_R, 
                                         constructs["R_gibson_ov_to_OLOH011", sequence]), .(enh_L, seq_L, enh_R, seq_R)]
seq[grepl("_B_", enh_L), refseq:= paste0(constructs["CASeq001_to_F_gibson_ov", sequence], seq_L, 
                                         constructs["rev_comp_TWIST_linker_R2", sequence],
                                         constructs["SCR1_ov", sequence], 
                                         constructs["TWIST_linker_F", sequence], seq_R, 
                                         constructs["R_gibson_ov_to_OLOH011", sequence]), .(enh_L, seq_L, enh_R, seq_R)]
seq[grepl("_C_", enh_L), refseq:= paste0(constructs["CASeq001_to_F_gibson_ov", sequence], seq_L, 
                                         constructs["rev_comp_TWIST_linker_R3", sequence],
                                         constructs["SCR1_ov", sequence], 
                                         constructs["TWIST_linker_F", sequence], seq_R, 
                                         constructs["R_gibson_ov_to_OLOH011", sequence]), .(enh_L, seq_L, enh_R, seq_R)]

# Create dat object
dat <- as.data.table(read_xlsx("metadata/sample_ids.xlsx"), colClasses = "character")
dat <- dat[order(Sample_ID)]
dat <- dat[, .(abfile_L= list.files(paste0("sanger_sequencing/sequencing_", Sequencing), paste0(Sample_ID, "_CASeq"), recursive = T, full.names = T),
               abfile_R= list.files(paste0("sanger_sequencing/sequencing_", Sequencing), paste0(Sample_ID, "_OLOH"), recursive = T, full.names = T)), c(colnames(dat))]
dat[seq, exp_seq:= i.refseq, on= c("ID_L", "ID_R")]

pdf("pdf/sanger_alignment_exp_sequences.pdf", height = 2, width = 10)
par(mfrow=c(1, 2), mar= c(2, 6, 2, 1), lwd= 10)
dat[, {obj <- my_sanger_align(abfile = abfile_L, refseq= exp_seq, features_v = constructs["SCR1_ov", sequence], features_names = "SPACER")
       my_sanger_plot(obj, labels = c(paste(ID_L, "vs", ID_R), paste(Sequencing, Sample_ID)))
       obj <- my_sanger_align(abfile = abfile_R, refseq= exp_seq, revcomp = T, features_v = constructs["SCR1_ov", sequence], features_names = "SPACER")
       my_sanger_plot(obj, labels = c(paste(ID_L, "vs", ID_R), paste(Sequencing, Sample_ID)))}, .(abfile_L, abfile_R, ID_L, ID_R, Sequencing, Sample_ID)]
dev.off()

# # Find closest sequences for the ones that dont align correctly
# bad <- data.table(Sample_ID= c("002", "022", "034", "218", "221", "222"),
#                   date= c("20200625", "20200625", "20200625", "20200703", "20200703", "20200703"))
# dat[bad, bad:= TRUE, on= c("Sample_ID", "date")]
# dat[bad=="TRUE", seq_L:= as.character(primarySeq(readsangerseq(abfile_L))), abfile_L]
# dat[bad=="TRUE", seq_R:= as.character(primarySeq(readsangerseq(abfile_R))), abfile_R]
# dat[bad=="TRUE" & nchar(seq_L>100), closest_L:= which.max(pairwiseAlignment(seq$refseq, seq_L, type = "global-local", scoreOnly= T)), seq_L]
# dat[bad=="TRUE" & nchar(seq_R>100), closest_R:= which.max(pairwiseAlignment(seq$refseq, as.character(reverseComplement(DNAString(seq_R))), type = "global-local", scoreOnly= T)), seq_R]
# dat[!is.na(closest_L), c("closest_ID_L_L", "closest_ID_R_L", "closest_enh_L_L", "closest_enh_R_L", "closest_seq_L"):= seq[closest_L, .(ID_L, ID_R, enh_L, enh_R, refseq)]]
# dat[!is.na(closest_R), c("closest_ID_L_R", "closest_ID_R_R", "closest_enh_L_R", "closest_enh_R_R", "closest_seq_R"):= seq[closest_R, .(ID_L, ID_R, enh_L, enh_R, refseq)]]
# saveRDS(dat, "Rdata/sanger_alignment.rds")
# 
# pdf("pdf/sanger_alignment_bad_sequences.pdf", height = 2, width = 20)
# par(mfrow=c(1, 2), mar= c(2, 6, 2, 1), lwd= 10)
# dat[!is.na(closest_seq_L) | !is.na(closest_seq_R), 
#   {
#      if(!is.na(closest_seq_L))
#      {
#        obj <- my_sanger_align(abfile = abfile_L, refseq= closest_seq_L, features_v = constructs["SCR1_ov", sequence], features_names = "SPACER")
#        my_sanger_plot(obj, labels = c(paste(closest_ID_L_L, "vs", closest_ID_R_L), paste(date, Sample_ID)))
#      }
#      if(!is.na(closest_seq_L))
#      {
#         obj <- my_sanger_align(abfile = abfile_R, refseq= closest_seq_R, revcomp = T, features_v = constructs["SCR1_ov", sequence], features_names = "SPACER")
#         my_sanger_plot(obj, labels = c(paste(closest_ID_L_R, "vs", closest_ID_R_R), paste(date, Sample_ID)))
#      }
#   }, .(abfile_L, abfile_R, closest_seq_L, closest_seq_R)]
# dev.off()

















