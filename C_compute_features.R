setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(GenomicRanges)
require(motifmatchr)
require(seqLogo)

lib <- readRDS("Rdata/uniq_library_final.rds")
lib <- lib[!detail=="ecoli"]

#-----------------------------------------------#
# 1- Gene assignment
#-----------------------------------------------#
# dm3 intervals lib file
enh <- GRanges(lib$coor, name= lib$ID)
enh <- enh[order(seqnames(enh), start(enh))]
# dm3 TSSs
tss <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene), 1, fix = "start")
tss$symbol <- mapIds(org.Dm.eg.db, key= tss$gene_id, column="SYMBOL", keytype="FLYBASE", multiVals="first")
tss <- tss[order(seqnames(tss), start(tss))]
seqlevels(tss, pruning.mode= "coarse") <- seqlevels(enh)
# Temp files
tmp1 <- tempfile(fileext = ".bed")
export(enh, tmp1)
tmp2 <- tempfile(fileext = ".bed")
export(tss, tmp2)
# All genes
all <- fread(cmd= closestBed(tmp1, tmp2))
all <- all[, .(V10= V10[1]), V1:V4]
# Add to lib
lib[all, ctss:= i.V10 , on= "ID==V4"]

#-----------------------------------------------#
# 2- dev STARR-Seq enrichment
#-----------------------------------------------#
# Compute gw STARR-Seq enrichment
STARR <- data.table(file= c("/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/input_DSCP_200bp_cut.bed",
                            "/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed", 
                            "/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed"), cdition= c("input", "screen", "screen"))
STARR <- STARR[, my_countReads(GRanges(lib$coor, name= lib$ID), file), c(colnames(STARR))]
# Filter out low input counts
check <- STARR[cdition=="input" & counts>75, name]
STARR <- STARR[name %in% check]
# Compute FE
STARR <- STARR[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
STARR <- dcast(STARR, name~cdition, value.var = "norm_counts")
STARR <- STARR[, .(name, DSCP200_log2FC= log(screen)-log2(input))]
# Add to lib
lib[STARR, DSCP200_log2FC:= i.DSCP200_log2FC, on= "ID==name"]

#-----------------------------------------------#
# 3- hk STARR-Seq enrichment
#-----------------------------------------------#
# Compute gw STARR-Seq enrichment
STARR <- data.table(file= c("/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/input_RPS12_200bp_cut.bed",
                            "/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/RpS12_200bp_gw_Rep1.UMI_cut.bed", 
                            "/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/RpS12_200bp_gw_Rep1.UMI_cut.bed"), cdition= c("input", "screen", "screen"))
STARR <- STARR[, my_countReads(GRanges(lib$coor, name= lib$ID), file), c(colnames(STARR))]
# Filter out low input counts
check <- STARR[cdition=="input" & counts>75, name]
STARR <- STARR[name %in% check]
# Compute FE
STARR <- STARR[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
STARR <- dcast(STARR, name~cdition, value.var = "norm_counts")
STARR <- STARR[, .(name, RPS200_log2FC= log(screen)-log2(input))]
# Add to lib
lib[STARR, RPS200_log2FC:= i.RPS200_log2FC, on= "ID==name"]

#-----------------------------------------------#
# 4- rep-STARR-Seq enrichment
#-----------------------------------------------#
# Compute gw sgl repSTARR-Seq enrichment
repSTARR <- data.table(file= c("/groups/stark/vloubiere/data/rep_STARRSeq_lorena/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_input.merged.uniq.bed",
                               "/groups/stark/vloubiere/data/rep_STARRSeq_lorena/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_1.UMI.bed", 
                               "/groups/stark/vloubiere/data/rep_STARRSeq_lorena/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_2.UMI.merged.bed"), cdition= c("input", "screen", "screen"))
repSTARR <- repSTARR[, my_countReads(GRanges(lib$coor, name= lib$ID), file), c(colnames(repSTARR))]
# Filter out low input counts
check <- repSTARR[cdition=="input" & counts>100, name]
repSTARR <- repSTARR[name %in% check]
# Compute FE
repSTARR <- repSTARR[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
repSTARR <- dcast(repSTARR, name~cdition, value.var = "norm_counts")
repSTARR <- repSTARR[, .(name, sgl_repSTARR_log2FC= log(screen)-log2(input))]
# Add to lib
lib[repSTARR, sgl_repSTARR_log2FC:= i.sgl_repSTARR_log2FC, on= "ID==name"]

#-----------------------------------------------#
# 5- Available data
#-----------------------------------------------#
# Compute ChIP-Seq ATAC-Seq enrichment
avail <- fread("/groups/stark/vloubiere/data/available_data/metadata/available_data_metadata.txt")
avail <- avail[sample!="input", .(file= uniq_bed_file, cdition= sample)]
avail <- avail[, my_countReads(GRanges(lib$coor, name= lib$ID), file, sorted = F), c(colnames(avail))]
# Filter out low input counts & compute enrichment
avail <- avail[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
avail[lib, control_group:= i.group=="control", on= "name==ID"]
avail[, control_norm_counts:= median(.SD[(control_group), norm_counts]), cdition]
avail[, log2FC:= log2(norm_counts)-log2(control_norm_counts)]
avail <- dcast(avail, name~paste0(cdition, "_log2FC"), value.var = "log2FC")
# Add to lib
lib <- lib[avail, , on= "ID==name"]

#-----------------------------------------------#
# 6- Motifs
#-----------------------------------------------#
som <- readRDS("Rdata/som_enriched_motifs.rds")
mot_cl <- as.data.table(t(as.data.table(som$codes)), keep.rownames = T)
colnames(mot_cl)[-1] <- paste0("motif_cl", seq(ncol(mot_cl)-1))
lib <- cbind(lib, mot_cl[, !"rn"])

saveRDS(lib, "Rdata/C_features.rds")
              
              
              