setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(seqinr)
require(gridExtra)
require(Biostrings)

#---------------------------------------#
# Retrieve all sequences
#---------------------------------------#
seq <- read.fasta("/groups/stark/vloubiere/genomes/Custom_peSTARRSeq_1/fasta/combinations_peSTARRSeq.fa", as.string = T)
constructs <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt")
setkeyv(constructs, "name")

# candidates
seq <- data.table(name= names(seq), seq= sapply(seq, function(x) as.character(x[1])))
candidate <- seq[!grep("^temp", name)]
temp_switch <- seq[grep("^temp", name)]

# plasmid
full_plasmid <- data.table(name= "plasmid", 
                           seq= constructs["DSCP_plasmid_full", 
                                           paste0(substr(sequence, start = nchar(sequence)-500, stop= nchar(sequence)), sequence, substr(sequence, start = 1, stop= 500))])
cloned_plasmid <- data.table(name= "cloned_plasmid",
                             seq= c(constructs[c("linker_R1_TWIST", "DSCP_downstream_cloning_site", "DSCP_upstream_cloning_site", "linker_F_TWIST"), paste0(sequence, collapse = "")],
                                    constructs[c("linker_R2_TWIST", "DSCP_downstream_cloning_site", "DSCP_upstream_cloning_site", "linker_F_TWIST"), paste0(sequence, collapse = "")],
                                    constructs[c("linker_R3_TWIST", "DSCP_downstream_cloning_site", "DSCP_upstream_cloning_site", "linker_F_TWIST"), paste0(sequence, collapse = "")]))

# constructs 
adapter <- data.table(name= c("left1_arm", "left2_arm", "left3_arm", "right_arm"), 
                      seq= c(constructs[c("linker_R1_TWIST", "rev_illumina_F", "CGCov_F"), paste0(sequence, collapse = "")],
                             constructs[c("linker_R2_TWIST", "rev_illumina_F", "CGCov_F"), paste0(sequence, collapse = "")],
                             constructs[c("linker_R3_TWIST", "rev_illumina_F", "CGCov_F"), paste0(sequence, collapse = "")],
                             constructs[c("CGCov_R", "rev_illumina_R", "linker_F_TWIST"), paste0(sequence, collapse = "")]))

# Spacers
spacer <- data.table(name= "SCR1_spacer", 
                     seq= constructs[c("rev_illumina_F", "SCR1_CGCov", "rev_illumina_R"), paste0(sequence, collapse="")])

dat <- rbindlist(list(candidate= candidate, temp_switch= temp_switch, adapter= adapter, spacer= spacer, full_plasmid= full_plasmid, cloned_plasmid= cloned_plasmid), idcol= T)
dat[, seq:= toupper(seq)]

#---------------------------------------#
# Restriction cuts
#---------------------------------------#
.e <- fread("/groups/stark/vloubiere/projects/z_miscellaneous/thermofisher_restriction_enzymes_table/thermofisher_restriction_enzymes_table.txt")
res <- dat[, .e[, .(name1, consensus_F, cutsite)], (dat)]
res[, match:= lengths(vmatchPattern(DNAString(consensus_F), DNAStringSet(seq), fixed= F))>0, .(name1, consensus_F)]$V1
pl <- res[, .(matches= length(which(match))), .(.id, name1, cutsite, consensus_F)]
pl <- dcast(pl, name1+cutsite+consensus_F~.id, value.var = "matches")
pl <- pl[order(candidate, spacer, adapter, temp_switch)]
.n <- unique(res[, .N, .(.id, name1, cutsite, consensus_F)][, .(.id, N)])
colnames(pl) <- .n[colnames(pl), , on=".id"][, paste0(.id, ifelse(is.na(N), "", paste0(" (", N, ")")))]

pdf("pdf/STARRSeq_design/restriction_sites_counts.pdf", width = 20, height = 55)
plot.new()
grid.table(pl)
mtext("Number cuts restriction enzyme")
dev.off()
