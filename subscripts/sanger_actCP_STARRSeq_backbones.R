setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
cdir <- "db/sanger_sequencing/220621_STARRSeq_backbones_actCPs/"
dir.create(cdir, showWarnings = F)

# files <- list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-06-22_PCW-IMP-232_0106/",
#                     "DSCPact",
#                     full.names = T)
# file.copy(files, cdir)

constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
dat <- data.table(file= list.files(cdir, 
                                   full.names = T, 
                                   recursive = T))
dat[, ID:= gsub("LOUB_(.*)_.*_.*$", "\\1", basename(file))]

# BglII/SbfI (SdaI) digest
dig_p <- vl_digest(constructs["DSCP_STARRSeq", sequence], c("BglII", "SdaI"))

dat[, seq:= paste0(dug_p[1], , dig_p[3])]


dat[file %like% "Fab7.*AatII", seq:= constructs["Fab7", sequence]]
dat[file %like% "Fab8.*AatII", seq:= constructs["Fab8", sequence]]
dat[file %like% "scs.*AatII", seq:= constructs["scs", sequence]]
dat[file %like% "scsprim.*AatII", seq:= constructs["scs'", sequence]]
dat[file %like% "gypsy.*AatII", seq:= constructs["gypsy", sequence]]
dat[, rev:= ifelse(grepl("CASeq001", file), F, T)]
dat[, seq:= paste0(dig_AatII[1], dig_AatII[2])]
# dat[, seq:= paste0(dig_AatII[1], seq, dig_AatII[2])]

dir.create("pdf/sanger_sequencing", showWarnings = F)
pdf("pdf/sanger_sequencing/220621_AatII.pdf")
dat[, vl_sanger_align(seq, 
                      abfiles = file, 
                      revcomp = rev, 
                      feat_sequences = constructs["sgl_LH", sequence], 
                      feat_cols = "black"), seq]
dev.off()
