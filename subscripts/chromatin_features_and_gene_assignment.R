setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

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
lib$SUHW_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_SuHw_rep1_uniq.bed"),
                                        peaks = lib[, .(seqnames, start, end)], 
                                        ext_peaks = 1000)$log2_enr

saveRDS(lib, "Rdata/300bp_enhancer_chromatin_features_and_gene_assignment.rds")
