setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(motifmatchr)
require(rtracklayer)
require(kohonen)
require(randomForest)

#------------------------------------------------------------------------------------------------------------------------------#
# Make enhancers unique and recompile features that were used for design 
#------------------------------------------------------------------------------------------------------------------------------#
vl8 <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist008_112019.rds"))
setnames(vl8, c("ID_vl", "enh_sequence"), c("ID", "enh_seq"))
vl12 <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))
lib <- rbindlist(list(T8= vl8, 
                      T12= vl12), 
                 fill= T)
lib <- lib[, .(seqnames, 
               start, 
               end, 
               strand,
               group, 
               detail, 
               enh_seq, 
               linker_ID, 
               ID)]
lib <- unique(lib)# Collapse the two libraries (some oligos are shared)
# Add values from BA screens
dev <- fread("../gw_STARRSeq_bernardo/db/peaks/DSCP_200bp_gw.UMI_cut_merged.peaks.txt")
dev[, c("seqnames", "start", "end"):= .(factor(V1), V2, V2+V3-1)]
lib$dev_log2FC_STARR200 <- dev[lib, ifelse(length(V7)>0, max(V7), as.numeric(NA)), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
hk <- fread("../gw_STARRSeq_bernardo/db/peaks/RpS12_200bp_gw.UMI_cut_merged.peaks.txt")
hk[, c("seqnames", "start", "end"):= .(factor(V1), V2, V2+V3-1)]
lib$hk_log2FC_STARR200 <- hk[lib, ifelse(length(V7)>0, max(V7), as.numeric(NA)), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
BA <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/BA_300bp_TWIST_STARRSeq.txt")
lib[BA, c("dev_log2FC_TWIST", "hk_log2FC_TWIST"):= .(i.dev_log2FoldChange, i.hk_log2FoldChange),
    on= c("seqnames", "start", "end", "strand")]
# Fix groups and details
lib[group=="Repressor", detail:= "Top_dip_LH", group]
lib[, group:= switch(as.character(group),
                     "ecdysone"= "inducible",
                     "heatshock"= "inducible",
                     "Repressor"= "repressor",
                     group), group]
lib[group %in% c("hk", "dev", "shared"), group:= fcase(dev_log2FC_TWIST>hk_log2FC_TWIST, "dev",
                                                       dev_log2FC_TWIST<=hk_log2FC_TWIST, "hk")]
lib[group %in% c("hk", "dev"), 
    detail:= cut(get(paste0(group, "_log2FC_TWIST")), 
                 c(-Inf,3,6,8,Inf),
                 include.lowest = T, 
                 labels = c("inactive", "weak", "medium", "strong")), group]
lib[group %in% c("hk", "dev") & abs(hk_log2FC_TWIST-dev_log2FC_TWIST)<1, group:= "shared"]
lib[, group:= factor(group,
                     levels= c("repressor",
                               "SUHW_peak",
                               "inducible",
                               "OSC",
                               "CP",
                               "control",
                               "DHS_peak",
                               "hk",
                               "shared",
                               "dev"))]
lib[, detail:= factor(detail,
                      levels= lib[, .N>0, detail][(V1), detail])]
# Add colors
lib[, col:= switch(as.character(group),
                   "hk" = "tomato",
                   "dev" = "#74C27A",
                   "control" = "lightgrey",
                   "OSC" = "black",
                   "inducible" = "gold",
                   "CP" = "deeppink2",
                   "DHS_peak" = "cyan3",
                   "shared" = "royalblue2",
                   "repressor" = "darkgrey",
                   "SUHW_peak" = "darkorchid3"), keyby= group]

#---------------------------------------------------------------#
# Gene assignment
#---------------------------------------------------------------#
# Closest promoter
gene <- import("/groups/stark/vloubiere/genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf")
seqlevelsStyle(gene) <- "UCSC"
tss <- as.data.table(gene)[type=="transcript", .(seqnames, start, end, strand, gene_id, symbol= gene_name)]
tss[, start:= ifelse(strand=="+", start, end)]
closestTss <- tss[lib, 
                  .(ID, 
                    closest_tss= symbol[which.min(abs(start-(i.start)))],
                    closest_tss_dist= min(abs(start-(i.start)))), 
                  .EACHI, 
                  on= "seqnames"][, .(ID, closest_tss, closest_tss_dist)]
# Multiple enhancer promoters
tss <- tss[symbol %in% c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo")]
closestSelTss <- tss[lib, 
                     .(ID, 
                       closest_sel_tss= symbol[which.min(abs(start-(i.start)))],
                       closest_sel_tss_dist= min(abs(start-(i.start)))),
                     .EACHI, on= "seqnames"][, .(ID, closest_sel_tss, closest_sel_tss_dist)]
closestTss <- merge(closestTss, 
                    closestSelTss)
lib <- merge(lib, closestTss, by= "ID")

#---------------------------------------------------------------#
# Chromatin Features
#---------------------------------------------------------------#
regions <- vl_resizeBed(bed= lib, 
                        center = "center", 
                        upstream = 1000, 
                        downstream = 1000, 
                        ignore.strand = T)
lib[, ATAC:= vl_bw_coverage(regions, "../available_data_dm3/db/bw/ATAC_merged.bw")]
lib[, SUHW:= vl_bw_coverage(regions, "../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw")]
lib[, H3K27Ac:= 
      vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K27ac_ChIP_Rep1.bw")+
      vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K27ac_ChIP_Rep2.bw")]
lib[, H3K4me1:= 
      vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K4me1_ChIP_Rep1.bw")+
      vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K4me1_ChIP_Rep2.bw")]
lib[, K4me3:= 
      vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE81795/tracks/S2_H3K4me3_Rep1.bw")+
      vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE81795/tracks/S2_H3K4me3_Rep2.bw")]
lib[, K27me3:= vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K27me3_ChIP.bw")]
lib[, PolII:= vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_PolII_ChIP.bw")]
lib[, GAF:= vl_bw_coverage(regions, "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE40646_Trl/tracks/S2_Trl_ChIP_Lis.bw")]

#--------------------------------------------------------------------#
# Deep STARR prediction
#--------------------------------------------------------------------#
# Script to run this: 
"/groups/stark/almeida/Projects/Vincent_pairedSTARRseq/results/DeepSTARR_predictions/Run_DeepSTARR.sh"
deep <- fread("/groups/stark/almeida/Projects/Vincent_pairedSTARRseq/results/DeepSTARR_predictions/vl_sequences.fa_predictions_Model_continuous_CNN_pool_120_7_60_3_60_3_pooling_2denses_dropout0.3_valloss2.73.txt")
setnames(deep, c("ID", "seq", "deep_dev", "deep_hk"))
lib <- merge(lib, deep[, .(ID, deep_dev, deep_hk)])

#------------------------------------------------------------------------------------------------------------------------------#
# TF motif counts
#------------------------------------------------------------------------------------------------------------------------------#
counts <- vl_motif_counts(sequences =  lib$enh_seq)
counts <- as.data.table(counts)
setnames(counts, function(x) paste0(x, "_motif"))
lib <- cbind(lib, counts)

cols <- grep("motif$", names(lib), value = T, invert = T)
setcolorder(lib, cols)
fwrite(lib,
       "Rdata/final_300bp_enhancer_features.txt", 
       sep= '\t',
       na = NA)
