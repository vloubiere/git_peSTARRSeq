#--------------------------------------------------------------------------####
# Compute enhancer features ####
lib <- as.data.table(readRDS("Rdata/uniq_300bp_enhancers.rds"))
# Closest promoter
gene <- import("/groups/stark/vloubiere/genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf")
seqlevelsStyle(gene) <- "UCSC"
tss <- as.data.table(gene)[type=="gene", .(seqnames, start, end, strand, gene_id, symbol= gene_name)]
tss[strand=="+", end:= start]
tss[strand=="-", start:= end]
tss[, ID:= paste0(seqnames, ":", start, "-", end, "__", symbol, "__", gene_id)]
lib$closest_tss <- tss[lib, ifelse(.N>0, 
                                   ID[which.min(abs(start-(i.start+144)))], 
                                   as.character(NA)), .EACHI, on= "seqnames"]$V1

# Multiple enhancer promoters
tss <- tss[c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo"), , on= "symbol"]
lib$closest_sel_tss <- tss[lib, ifelse(.N>0 & min(abs(start-(i.start+144)))<100000, 
                                       ID[which.min(abs(start-(i.start+144)))], 
                                       as.character(NA)), .EACHI, on= "seqnames"]$V1

# STARR-Seq/ChIP-Seq enrichment
lib$sgl_STARR_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_1.UMI.bed", 
                                                          "/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_2.UMI.merged.bed"),
                                             Input_bed = c("/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_input.merged.uniq.bed"), 
                                             peaks = lib[, .(seqnames, start, end)], 
                                             ext_peaks = 1000)$log2_enr
lib$ATAC_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep2_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep3_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep4_uniq.bed"),
                                        peaks = lib[, .(seqnames, start, end)], 
                                        ext_peaks = 1000)$log2_enr
lib$GAF_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep1_uniq.bed",
                                                    "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep2_uniq.bed"),
                                       Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep2_uniq.bed"),
                                       peaks = lib[, .(seqnames, start, end)], 
                                       ext_peaks = 1000)$log2_enr
lib$H3K27ac_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep1_uniq.bed",
                                                        "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep2_uniq.bed"),
                                           Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"),
                                           peaks = lib[, .(seqnames, start, end)], 
                                           ext_peaks = 1000)$log2_enr
lib$H3K4me1_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep1_uniq.bed",
                                                        "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep2_uniq.bed"),
                                           Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"),
                                           peaks = lib[, .(seqnames, start, end)], 
                                           ext_peaks = 1000)$log2_enr
lib$H3K4me3_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep1_uniq.bed",
                                                        "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep2_uniq.bed"),
                                           Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep2_uniq.bed"),
                                           peaks = lib[, .(seqnames, start, end)], 
                                           ext_peaks = 1000)$log2_enr
lib$Pol2_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_RNAPolII_rep1_uniq.bed"),
                                        Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed"),
                                        peaks = lib[, .(seqnames, start, end)], 
                                        ext_peaks = 1000)$log2_enr
lib$H3K27me3_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27me3_rep1_uniq.bed"),
                                            Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed"),
                                            peaks = lib[, .(seqnames, start, end)], 
                                            ext_peaks = 1000)$log2_enr


#---------------------------------#
# Compute counts for most informative motifs
#---------------------------------#
som <- readRDS("Rdata/som_clustering_motifs_300bp_enhancers.rds")
setorderv(som$info, "dist")
sel <- som$info[!is.na(BA_cluster), .SD[1, .(motif, Dmel, BA_cluster)], cl]
sel_idx <- name(TF_clusters_PWMs$All_pwms_log_odds) %in% sel$motif
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel_idx], 
                   GRanges(lib[!grepl("Ecoli", ID), .(seqnames, start, end)]), 
                   genome= "dm3", 
                   p.cutoff= 5e-4, 
                   bg="even", 
                   out= "scores")
counts <- as.matrix(motifCounts(hit))
colnames(counts) <- paste0("motif__", name(TF_clusters_PWMs$All_pwms_log_odds)[sel_idx])
rownames(counts) <- lib[!grepl("Ecoli", ID), ID]
counts <- as.data.table(counts, keep.rownames = T)

lib <- merge(lib, 
             counts, 
             by.x= "ID",
             by.y= "rn", 
             all.x= T)


# SAVE
setcolorder(lib, c("ID", 
                   "group", 
                   "detail", 
                   "vl", 
                   "linker_ID", 
                   "col", 
                   "seqnames", 
                   "start", 
                   "end", 
                   "strand", 
                   "closest_tss", 
                   "closest_sel_tss",
                   "dev_log2FoldChange", 
                   "hk_log2FoldChange"))
saveRDS(lib, "Rdata/final_300bp_enhancer_features.rds")










