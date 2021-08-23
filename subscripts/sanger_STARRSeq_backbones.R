setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir.create("db/sanger_sequencing/220621_STARRSeq_backbones_MP/", showWarnings = F)

# files <- list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-06-22_PCW-IMP-232_0106/",
#                     "_pGL3DSCPII_|_pGL3RpS12_",
#                     full.names = T)
# file.copy(files, "db/sanger_sequencing/220621_STARRSeq_backbones_MP/")

# source("/groups/stark/vloubiere/exp_data/update_files.R")

constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")

#--------------------------------#
dat <- data.table(file= list.files("db/sanger_sequencing/220621_STARRSeq_backbones_MP/", 
                                   full.names = T, 
                                   recursive = T))
dat[, ID:= gsub("LOUB_(.*)_.*_.*$", "\\1", basename(file))]
DSCP_seq <- constructs["DSCP_STARRSeq", sequence]
dat[ID %like% "DSCPII", seq:= paste0(substr(DSCP_seq, 4500, nchar(DSCP_seq)),
                                     substr(DSCP_seq, 1, 1500))]
RPS_seq <- constructs["RPS12_STARRSeq", sequence]
dat[ID %like% "RpS12", seq:= paste0(substr(RPS_seq, 4500, nchar(RPS_seq)),
                                    substr(RPS_seq, 1, 1500))]
dat[, rev:= ifelse(grepl("CASeq001", file), F, T)]

#--------------------------------#
dir.create("pdf/sanger_sequencing", showWarnings = F)

pdf("pdf/sanger_sequencing/220621_STARRSeq_backbones_MP.pdf")
dat[, vl_sanger_align(seq, 
                      abfiles = file, 
                      revcomp = rev, 
                      feat_sequences = constructs["sgl_LH", sequence], 
                      feat_cols = "black"), seq]
dev.off()
