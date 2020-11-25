setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(seqinr)
require(gridExtra)

#---------------------------------------#
# Retrieve all sequences
#---------------------------------------#
# candidates
seq <- read.fasta("/groups/stark/vloubiere/genomes/Custom_peSTARRSeq_1/fasta/combinations_peSTARRSeq.fa", as.string = T)
seq <- data.table(name= names(seq), seq= sapply(seq, function(x) as.character(x[1])))
candidate <- seq[!grep("^temp", name)]
temp_switch <- seq[grep("^temp", name)]
# Adapters 
constructs <- fread("Rdata/constructs_sequences/constructs_sequences.txt")
setkeyv(constructs, "name")
adapter <- data.table(name= c("left_arm", "right1_arm", "right2_arm", "right3_arm"), 
                      seq= c(constructs[c("DSCP_upstream_from_CASEq044_to_GIBSONov", "illumina_F", "linker_F_TWIST"), 
                                        paste0(sequence, collapse = "")],
                             constructs[c("linker_R1_TWIST", "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), 
                                        paste0(sequence, collapse = "")],
                             constructs[c("linker_R2_TWIST", "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), 
                                        paste0(sequence, collapse = "")],
                             constructs[c("linker_R3_TWIST", "illumina_R", "DSCP_downstream_from_GIBSONov_to_CASeq003"), 
                                        paste0(sequence, collapse = "")]))
# Spacers
spacer <- data.table(name= paste0("SCR1_spacer_", 1:3), 
                     seq= c(constructs[c("linker_R1_TWIST", "SCR1_CGCov", "linker_F_TWIST"), paste0(sequence, collapse="")],
                            constructs[c("linker_R2_TWIST", "SCR1_CGCov", "linker_F_TWIST"), paste0(sequence, collapse="")],
                            constructs[c("linker_R3_TWIST", "SCR1_CGCov", "linker_F_TWIST"), paste0(sequence, collapse="")]))
dat <- rbindlist(list(candidate= candidate, temp_switch= temp_switch, adapter= adapter, spacer= spacer), idcol= T)
dat[, seq:= toupper(seq)]

#---------------------------------------#
# Restriction cuts
#---------------------------------------#
enz <- data.table(enzyme= c("HindIII", "KasI", "MreI"), 
                  site= c("AAGCTT|TTCGAA", "GGCGCC|CCGCGG", "CGCCGGCG|GCGGCCGC"))
dat <- dat[, (enz), (dat)]
dat[, cuts:= grepl(site, seq), .(site, seq)]

res <- dat[, .(N_cuts= paste0(length(which(cuts)), "/", .N)), .(.id, enzyme, site)]
res <- dcast(res, .id~enzyme+site, value.var = "N_cuts")

pdf("pdf/intron_STARRSeq/restriction_sites_counts.pdf", width = 10, height = 5)
plot.new()
grid.table(res)
mtext("Number cuts restriction enzyme")
dev.off()






