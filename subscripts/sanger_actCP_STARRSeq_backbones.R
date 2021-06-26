setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
cdir <- "db/sanger_sequencing/220621_STARRSeq_backbones_actCPs/"
dir.create(cdir, showWarnings = F)

# files <- list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-06-22_PCW-IMP-232_0106/",
#                     "DSCPact",
#                     full.names = T)
# files <- list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-06-24_PCW-IMP-232_0117/",
#                     "DSCPact",
#                     full.names = T)
# file.copy(files, cdir)

# source("/groups/stark/vloubiere/exp_data/update_files.R")

constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
dat <- data.table(file= list.files(cdir, 
                                   "CASeq001|CASeq002",
                                   full.names = T, 
                                   recursive = T))
dat[, ID:= gsub("LOUB_(.*)_.*_.*$", "\\1", basename(file))]

# BglII/SbfI (SdaI) digest
dig_p <- vl_digest(constructs["DSCP_STARRSeq", sequence], 
                   c("BglII", "SdaI"), 
                   keepsite = T)
dat[grepl("DSCPact1", ID), seq:= paste0(substr(dig_p[3], nchar(dig_p[3])-1500, nchar(dig_p[3])),
                                        dig_p[1], 
                                        constructs["actCP_1", sequence], 
                                        substr(dig_p[3], 1, 1500))]
dat[grepl("DSCPact2", ID), seq:= paste0(substr(dig_p[3], nchar(dig_p[3])-1500, nchar(dig_p[3])),
                                        dig_p[1], 
                                        constructs["actCP_2", sequence], 
                                        substr(dig_p[3], 1, 1500))]
# dat[, seq:= paste0(dig_p[1],dig_p[3])]
dat[, rev:= ifelse(grepl("CASeq001", file), F, T), file]

pdf("pdf/sanger_sequencing/actCP_STARRSeq_backbones.pdf")
dat[, vl_sanger_align(seq, 
                      abfiles = file, 
                      revcomp = rev, 
                      feat_sequences = constructs["actCP_2", sequence], 
                      feat_cols = "black"), seq]
dev.off()
