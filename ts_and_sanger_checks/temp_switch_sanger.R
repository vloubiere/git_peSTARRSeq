setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

seq <- fread("/groups/stark/vloubiere/exp_data/constructs_sequences.txt")
setkeyv(seq, "name")

dat <- data.table(name= c("HAM1", "SCR2", "SUP1"),
                  sequence= list(seq[c("DSCP_upstream_from_CASEq044_to_GIBSONov", "illumina_F", 
                                       "linker_F_TWIST", "HAM1", "linker_R1_TWIST", 
                                       "SCR1_CGCov", 
                                       "linker_F_TWIST", "HAM1", "linker_R1_TWIST", 
                                       "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), sequence],
                                 seq[c("DSCP_upstream_from_CASEq044_to_GIBSONov", "illumina_F", 
                                       "linker_F_TWIST", "SCR2", "linker_R1_TWIST", 
                                       "SCR1_CGCov", 
                                       "linker_F_TWIST", "SCR2", "linker_R1_TWIST", 
                                       "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), sequence],
                                 seq[c("DSCP_upstream_from_CASEq044_to_GIBSONov", "illumina_F", 
                                       "linker_F_TWIST", "SUP1", "linker_R1_TWIST", 
                                       "SCR1_CGCov", 
                                       "linker_F_TWIST", "SUP1", "linker_R1_TWIST", 
                                       "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), sequence]))

dat[, file:= .(list(list.files("db/sanger_sequencing/libvl002_singcol_bulk_tempswitch", name, full.names = T))), name]
dat[, revcomp:= lapply(file, function(x) ifelse(grepl("CASeq044", x), "F", "T"))]

pdf("pdf/ts_and_sanger_checks/template_switching.pdf", height = 6, width = 9)
par(mfrow= c(3,1))
dat[, my_sanger_align(refseq = unlist(sequence), abfiles_vec =  unlist(file), main= name[1], 
                      revcomp = unlist(revcomp), refseq_names = c("pl", "", "", name, "", "SPACER", "", name, "", "", "pl")), name]
dev.off()

