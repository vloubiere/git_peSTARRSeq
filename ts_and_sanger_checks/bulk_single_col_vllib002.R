setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

seq <- fread("Rdata/constructs_sequences/constructs_sequences.txt")
setkeyv(seq, "name")

dat <- data.table(file= list.files("db/sanger_sequencing/libvl002_singcol_bulk_tempswitch", "bulk|singcol", full.names = T))
dat[, name:= tstrsplit(basename(file), "_", keep= 2)]
dat <- dat[, .(file= list(file), refseq= list(seq[c("DSCP_upstream_from_CASEq044_to_GIBSONov", "illumina_F", 
                                                    "linker_F_TWIST", "SCR2", "linker_R1_TWIST", 
                                                    "SCR1_CGCov", 
                                                    "linker_F_TWIST", "SCR2", "linker_R1_TWIST", 
                                                    "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), sequence]),
               refseq_names= list(c("pl", "", "", "SCR2", "", "SPACER", "", "SCR2", "", "", "pl"))), name]
dat[, revcomp:= lapply(file, function(x) ifelse(grepl("CASeq044", x), "F", "T"))]

pdf("pdf/ts_and_sanger_checks/bulk_sing_col_vllib002.pdf", width = 25, height = 15)
par(mfrow= c(7,3))
dat[, my_sanger_align(refseq = unlist(refseq), abfiles_vec =  unlist(file), main= name[1], 
                      revcomp = unlist(revcomp), refseq_names = unlist(refseq_names)), name]
dev.off()