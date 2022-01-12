setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(motifmatchr)
require(rtracklayer)

#---------------------------------------------------------------#
# Clean object containing available data
#---------------------------------------------------------------#
# Merge BA and my TWIST into a single non-redundant clean object
vl8 <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
names(vl8)[names(vl8)=="ID_vl"] <- "ID"
vl12 <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))
lib <- rbind(vl8, vl12, fill= T)
lib <- lib[, .(seqnames, 
               start, 
               end, 
               strand, 
               group, 
               detail, 
               oligo_full_sequence, 
               linker_ID, 
               ID)]
# Collapse the two libraries (some oligos are shared)
lib <- unique(lib)
# Add values from BA screen
BA <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
lib[, c("dev_log2FoldChange", "hk_log2FoldChange"):= 
      BA[.BY, .(dev_log2FoldChange, hk_log2FoldChange), 
         on= c("seqnames", "start", "end", "strand")], .(seqnames, start, end, strand)]
# Uniq groups
lib[group=="shared", group:= ifelse(dev_log2FoldChange>hk_log2FoldChange, "dev", "hk")]
lib[group=="ecdysone", c("group", "detail"):= .("inducible", "ecdysone")]
lib[group=="heatshock", c("group", "detail"):= .("inducible", "heatshock")]
lib[group=="Repressor", group:= "repressor"]
lib[, group:= factor(group, 
                     c("repressor", "SUHW_peak", "inducible", "OSC", "CP", "control", "DHS_peak", "dev", "hk"))]
# Uniq details
lib[group== "dev", 
    detail:= cut(dev_log2FoldChange, c(-Inf,4,6,8,Inf), labels = c("inactive", "weak", "medium", "strong"))]
lib[group== "hk", 
    detail:= cut(hk_log2FoldChange, 
                 quantile(hk_log2FoldChange, seq(0, 1, length.out= 5)), 
                 include.lowest = T, 
                 labels = c("inactive", "weak", "medium", "strong"))]
lib[group=="control" & detail=="Ecoli", detail:= "ecoli"]
lib[, detail:= factor(detail, 
                      as.data.table(table(lib$detail))[N>0, V1])]
# Add colors
class_Cc <- data.table(group= c("hk", "dev", "OSC", "inducible", "control", "CP", "DHS_peak", "repressor", "SUHW_peak"),
                       col= c("tomato", "#74C27A", "black", "gold", "lightgrey", "cyan3", "royalblue2", "deeppink2", "darkorchid3"))
lib <- lib[class_Cc, , on= "group"]

#---------------------------------------------------------------#
# Gene assignment
#---------------------------------------------------------------#
# Closest promoter
gene <- import("/groups/stark/vloubiere/genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf")
seqlevelsStyle(gene) <- "UCSC"
tss <- as.data.table(gene)[type=="transcript", .(seqnames, start, end, strand, gene_id, symbol= gene_name)]
tss[, start:= ifelse(strand=="+", start, end)]
lib$closest_tss <- tss[lib, ifelse(.N>0, 
                                   symbol[which.min(abs(start-(i.start)))], 
                                   as.character(NA)), .EACHI, on= "seqnames"]$V1

# Multiple enhancer promoters
tss <- tss[c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo"), , on= "symbol"]
lib$closest_sel_tss <- tss[lib, ifelse(.N>0 & min(abs(start-i.start))<100000, 
                                       symbol[which.min(abs(start-i.start))], 
                                       as.character(NA)), .EACHI, on= "seqnames"]$V1

#---------------------------------------------------------------#
# DNA motifs clustering
#---------------------------------------------------------------#
# Count hits low cutoff
hits <- vl_motif_counts(sequences =  lib$oligo_full_sequence)
# Select motifs enriched in at least one of the lib groups
enr <- vl_motif_cl_enrich(hits, 
                          cl_IDs = lib$group, 
                          plot = F, 
                          N_top = 10)
motifs <- unique(enr[padj<0.00001 & log2OR>1, motif])
#---------------------------------------------------------------#
# Compute motifs'clusters
#---------------------------------------------------------------#
cl <- hits[, ..motifs]
# Clip outliers and log
cl <- cl[, lapply(.SD, function(x) {
  lim <- quantile(x, c(0.05, 0.95))
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  log2(x+1)
})]
# Cluster and add to lib objectsave heatmap clustering
final_cl <- vl_heatmap(as.matrix(cl), 
                       cutree_rows = 12,
                       clustering_distance_cols = "spearman",
                       plot = F)
pdf("pdf/motifs_clustering/heatmap_motif_clustering_300bp_uniq_enhancers.pdf", 
    height = 10)
plot(final_cl,
     breaks = c(0, 3), 
     col = c("blue", "yellow"), 
     cutree_rows= 12,
     cutree_cols = 12)
dev.off()
lib[, motif_cl:= final_cl$result_DT[, unique(rn_cl), keyby= as.numeric(rn)]$V1]
top_motifs <- unique(enr[padj<1e-10 & log2OR>1, motif])
top_motifs <- hits[, ..top_motifs][, lapply(.SD, function(x) log2(x+1))]
names(top_motifs) <- paste0("top_motif__", names(top_motifs))
lib <- cbind(lib, top_motifs)

#---------------------------------------------------------------#
# Chromatin features
#---------------------------------------------------------------#
# ChIP-Seq enrichment
if(!file.exists("db/bed/GSE119708_ATAC_merged.bed"))
{
  ATAC <- list.files("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/",
                      "ATAC_rep.*.bed$", 
                     full.names = T)
  system(paste(c("cat", ATAC, ">", "db/bed/GSE119708_ATAC_merged.bed"), collapse = " "))
  vl_exportBed(vl_shuffleBed("db/bed/GSE119708_ATAC_merged.bed"), "db/bed/GSE119708_ATAC_shuffled.bed")
}
lib$ATAC_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                            center = "center", 
                                                            upstream = 500, 
                                                            downstream = 500, 
                                                            ignore.strand = T),
                                     ChIP_bed = "db/bed/GSE119708_ATAC_merged.bed",
                                     Input_bed = "db/bed/GSE119708_ATAC_shuffled.bed")$signalValue)

lib$H3K27ac_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                               center = "center", 
                                                               upstream = 500, 
                                                               downstream = 500, 
                                                               ignore.strand = T),
                                        ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep2_uniq.bed"),
                                        Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"))$signalValue)

lib$H3K4me1_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                               center = "center", 
                                                               upstream = 500, 
                                                               downstream = 500, 
                                                               ignore.strand = T),
                                        ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep2_uniq.bed"),
                                        Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"))$signalValue)

lib$H3K4me3_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                               center = "center", 
                                                               upstream = 500, 
                                                               downstream = 500, 
                                                               ignore.strand = T),
                                        ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep2_uniq.bed"),
                                        Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep1_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep2_uniq.bed"))$signalValue)

lib$H3K27me3_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                center = "center", 
                                                                upstream = 500, 
                                                                downstream = 500, 
                                                                ignore.strand = T),
                                         ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27me3_rep1_uniq.bed"),
                                         Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed"))$signalValue)

lib$Pol2_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                center = "center", 
                                                                upstream = 500, 
                                                                downstream = 500, 
                                                                ignore.strand = T),
                                     ChIP_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_RNAPolII_rep1_uniq.bed",
                                     Input_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed")$signalValue)

lib$GAF_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                            center = "center", 
                                                            upstream = 500, 
                                                            downstream = 500, 
                                                            ignore.strand = T),
                                    ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep1_uniq.bed",
                                                 "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep2_uniq.bed"),
                                    Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep1_uniq.bed",
                                                  "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep2_uniq.bed"))$signalValue)

lib$SUHW_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                           center = "center", 
                                                           upstream = 500, 
                                                           downstream = 500, 
                                                           ignore.strand = T),
                                    ChIP_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_SuHw_rep1_uniq.bed",
                                    Input_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_input_rep1_uniq.bed")$signalValue)

#------------------------------------------------------------------------------------------------#
# SAVE
#------------------------------------------------------------------------------------------------#
saveRDS(lib, "Rdata/final_300bp_enhancer_features.rds")
