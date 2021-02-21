setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(GenomicRanges)
require(seqinr)
require(motifmatchr)
require(PWMEnrich)
require(TFBSTools)
require(seqLogo)

lib <- readRDS("Rdata/library/uniq_library_final.rds")
lib <- lib[detail!="ecoli"]

#-----------------------------------------------#
# 1- Gene assignment
#-----------------------------------------------#
# dm3 intervals lib file
enh <- GRanges(lib$coor, name= lib$ID)
enh <- enh[order(as.character(seqnames(enh)), start(enh))]
# dm3 TSSs
tss <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene), 1, fix = "start")
tss$symbol <- mapIds(org.Dm.eg.db, key= tss$gene_id, column="SYMBOL", keytype="FLYBASE", multiVals="first")
tss <- tss[order(as.character(seqnames(tss)), start(tss))]
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
STARR <- data.table(file= c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/input_DSCP_200bp_cut.bed",
                            "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed", 
                            "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed"), 
                    cdition= c("input", "screen", "screen"))
STARR <- STARR[, my_countReads(GRanges(lib$coor, name= lib$ID), file), c(colnames(STARR))]
# Filter out low input counts
check <- STARR[cdition=="input" & counts>75, name]
STARR <- STARR[name %in% check]
# Compute FE
STARR <- STARR[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
STARR <- data.table::dcast(STARR, name~cdition, value.var = "norm_counts")
STARR <- STARR[, .(name, DSCP200_log2FoldChange= log2(screen)-log2(input))]
# Add to lib
lib[STARR, DSCP200_log2FoldChange:= i.DSCP200_log2FoldChange, on= "ID==name"]

#-----------------------------------------------#
# 3- hk STARR-Seq enrichment
#-----------------------------------------------#
# Compute gw STARR-Seq enrichment
STARR <- data.table(file= c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/input_RPS12_200bp_cut.bed",
                            "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/RpS12_200bp_gw_Rep1.UMI_cut.bed", 
                            "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/RpS12_200bp_gw_Rep1.UMI_cut.bed"), 
                    cdition= c("input", "screen", "screen"))
STARR <- STARR[, my_countReads(GRanges(lib$coor, name= lib$ID), file), c(colnames(STARR))]
# Filter out low input counts
check <- STARR[cdition=="input" & counts>75, name]
STARR <- STARR[name %in% check]
# Compute FE
STARR <- STARR[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
STARR <- data.table::dcast(STARR, name~cdition, value.var = "norm_counts")
STARR <- STARR[, .(name, RPS200_log2FoldChange= log2(screen)-log2(input))]
# Add to lib
lib[STARR, RPS200_log2FoldChange:= i.RPS200_log2FoldChange, on= "ID==name"]

#-----------------------------------------------#
# 4- rep-STARR-Seq enrichment
#-----------------------------------------------#
# Compute gw sgl repSTARR-Seq enrichment
repSTARR <- data.table(file= c("/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_input.merged.uniq.bed",
                               "/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_1.UMI.bed", 
                               "/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_2.UMI.merged.bed"), 
                       cdition= c("input", "screen", "screen"))
repSTARR <- repSTARR[, my_countReads(GRanges(lib$coor, name= lib$ID), file), c(colnames(repSTARR))]
# Filter out low input counts
check <- repSTARR[cdition=="input" & counts>100, name]
repSTARR <- repSTARR[name %in% check]
# Compute FE
repSTARR <- repSTARR[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
repSTARR <- data.table::dcast(repSTARR, name~cdition, value.var = "norm_counts")
repSTARR <- repSTARR[, .(name, sgl_repSTARR_log2FoldChange= log2(screen)-log2(input))]
# Add to lib
lib[repSTARR, sgl_repSTARR_log2FoldChange:= i.sgl_repSTARR_log2FoldChange, on= "ID==name"]

#-----------------------------------------------#
# 5- Available data
#-----------------------------------------------#
# Compute ChIP-Seq ATAC-Seq enrichment
avail <- fread("/groups/stark/vloubiere/projects/available_data/db/metadata/available_data_metadata.txt")
avail[, file:= list.files("/groups/stark/vloubiere/projects/available_data/db/bed/", paste0(uniq_name, ".*.bed"), full.names = T), uniq_name]
avail <- avail[sample!="input", .(file, cdition= sample)]
avail <- avail[, my_countReads(GRanges(lib$coor, name= lib$ID), file, sorted = F), c(colnames(avail))]
# Filter out low input counts & compute enrichment
avail <- avail[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), c("name", "cdition")]
avail[lib, control_group:= i.group=="control", on= "name==ID"]
avail[, control_norm_counts:= median(.SD[(control_group), norm_counts]), cdition]
avail[, log2FC:= log2(norm_counts)-log2(control_norm_counts)]
avail <- data.table::dcast(avail, name~paste0(cdition, "_log2FC"), value.var = "log2FC")
# Add to lib
lib <- lib[avail, , on= "ID==name"]

#-----------------------------------------------#
# 6- Motifs
#-----------------------------------------------#
# Import som clustering and identify the most informative/representative motifs
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
som <- readRDS("Rdata/motifs/som_enriched_motifs.rds")
sel <- match(som$info$best_match, name(TF_clusters_PWMs$All_pwms_log_odds))

# counts
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], GRanges(lib$coor), genome= "dm3", p.cutoff= 5e-4, bg="even", out= "scores")
counts <- as.matrix(motifCounts(hit))
colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds[sel])
rownames(counts) <- lib$ID
counts <- as.data.table(counts, keep.rownames = T)
colnames(counts)[-1] <- paste0("motif__", colnames(counts)[-1])

# add to lib and add ecoli sequences
lib <- cbind(lib, counts[, !"rn"])
lib <- rbind(lib, 
             data.table(ID= grep("Ecoli", readRDS("Rdata/library/vl_library_112019.rds")$ID_vl, value= T), 
                        group= "control", detail= "ecoli", col= "lightgrey"), fill= T)

#-------------------------------------------#
# 7- Add Bernardo's clusters
#-------------------------------------------#
source("git_peSTARRSeq/peSTARRSeq_data_processing/_som_motifs_clustering_enhancer_groups.R")
lib <- addClustersToLib(lib)

# SAVE
saveRDS(lib, "Rdata/library/lib_features.rds")

              