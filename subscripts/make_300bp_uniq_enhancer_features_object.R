setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(motifmatchr)
require(rtracklayer)
require(kohonen)
require(randomForest)

#------------------------------------------------------------------------------------------------------------------------------#
# Make enhancers unique and recompile features that were used for design 
#------------------------------------------------------------------------------------------------------------------------------#
# Merge the two TWIST libraries into a single non-redundant object
if(!file.exists("Rdata/uniq_enh_feat/lib.rds"))
{
  vl8 <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist008_112019.rds"))
  setnames(vl8, "ID_vl", "ID")
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
                 oligo_full_sequence, 
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
  saveRDS(lib, 
          "Rdata/uniq_enh_feat/lib.rds")
}else
  lib <- readRDS( "Rdata/uniq_enh_feat/lib.rds")

#------------------------------------------------------------------------------------------------------------------------------#
# TF motif 
#------------------------------------------------------------------------------------------------------------------------------#
# Count hits low cutoff
if(!file.exists("Rdata/uniq_enh_feat/hits.rds"))
{
  set.seed(1)
  rdm <- vl_control_regions_BSgenome(BSgenome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3, 
                                     n = 5000, 
                                     width = 299)
  rdm <-BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3, 
                         rdm$seqnames, 
                         rdm$start, 
                         rdm$end, 
                         as.character= T)
  hits <- vl_motif_counts(sequences =  c(lib$oligo_full_sequence, rdm))
  setnames(hits, function(x) vl_Dmel_motifs_DB_full[x, uniqName_noSpecialChar, on= "motif"])
  saveRDS(hits, 
          "Rdata/uniq_enh_feat/hits.rds")
}else
  hits <- readRDS("Rdata/uniq_enh_feat/hits.rds")

# Select motifs enriched in at least one of the lib groups
if(!file.exists("Rdata/uniq_enh_feat/motifEnrichOverRdm.rds"))
{
  enr <- vl_motif_cl_enrich(as.matrix(hits),
                            cl_IDs = c(rep("lib", nrow(lib)),
                                       rep("rdm", length(rdm))))
  saveRDS(enr,
          "Rdata/uniq_enh_feat/motifEnrichOverRdm.rds")
  
}else
  enr <- readRDS("Rdata/uniq_enh_feat/motifEnrichOverRdm.rds")

if(!file.exists("Rdata/uniq_enh_feat/motifCounts.rds"))
{
  motifs <- as.character(unique(enr[cl=="lib" & padj<0.05 & log2OR>0, motif]))
  motifCounts <- hits[1:nrow(lib), ..motifs]
  motifCounts[, ID:= lib$ID]
  setcolorder(motifCounts, "ID")
  saveRDS(motifCounts,
          "Rdata/uniq_enh_feat/motifCounts.rds")
}else
  motifCounts <- readRDS("Rdata/uniq_enh_feat/motifCounts.rds")

#---------------------------------------------------------------#
# Compute motifs'clusters
#---------------------------------------------------------------#
if(!file.exists("Rdata/uniq_enh_feat/somMotifsClustering.rds"))
{
  cap <- 0.95
  mat <- apply(as.matrix(motifCounts, 1), 2, function(x) {
    x <- ifelse(x>quantile(x, cap), quantile(x, cap), x)
    return(scale(x))
  })
  mat <- mat[,!apply(mat, 2, anyNA)]
  grid <- somgrid(xdim = 10,
                  ydim = 10,
                  topo = "hexagonal", 
                  toroidal = T)
  set.seed(1)
  som <- supersom(mat, 
                  grid = grid)
  saveRDS(som, "Rdata/uniq_enh_feat/somMotifsClustering.rds")
}else
  som <- readRDS("Rdata/uniq_enh_feat/somMotifsClustering.rds")

if(!file.exists("Rdata/uniq_enh_feat/somMotifsClustering_kmeanNodes.rds"))
{
  set.seed(1)
  km <- kmeans(som$codes[[1]], 9)
  saveRDS(km,
          "Rdata/uniq_enh_feat/somMotifsClustering_kmeanNodes.rds")
}else
  km <- readRDS("Rdata/uniq_enh_feat/somMotifsClustering_kmeanNodes.rds")
km_cl <- km$cluster[som$unit.classif]

# Check result
heat <- dcast(lib[, .(group, cl= paste0("cl_", km_cl))],
              cl~group,
              fun.aggregate = length)
vl_heatmap(scale(as.matrix(heat, 1)),
           breaks = c(-2,0,2))

# Top motifs based on clustering
if(!file.exists("Rdata/uniq_enh_feat/motifEnrichKmClusters.rds"))
{
  enr_cl <- vl_motif_cl_enrich(counts_matrix = as.matrix(hits),
                               cl_IDs = c(km_cl, rep("rdm", length(rdm))), 
                               control_cl = "rdm")
  setorderv(enr_cl, c("cl", "log2OR", "padj"), order = c(1,-1,1))
  enr_cl[is.finite(log2OR), rank:= seq(.N), cl]
  saveRDS(enr_cl,
          "Rdata/uniq_enh_feat/motifEnrichKmClusters.rds")
}else
  enr_cl <- readRDS("Rdata/uniq_enh_feat/motifEnrichKmClusters.rds")

if(!file.exists("Rdata/uniq_enh_feat/topMotifCountsKmClusters.rds"))
{
  top_motifs <- as.character(unique(enr_cl[rank<25, motif]))
  topMotifCounts <- hits[1:nrow(lib), ..top_motifs]
  saveRDS(topMotifCounts,
          "Rdata/uniq_enh_feat/topMotifCountsKmClusters.rds")
}else
  topMotifCounts <- readRDS("Rdata/uniq_enh_feat/topMotifCountsKmClusters.rds")
  
# SOM clustering in order to get codes(~100Go of RAM!!). Useful to avoid modelling overfitting
if(!file.exists("Rdata/uniq_enh_feat/som_codes_enh_motifs_codesOnly_n16.rds"))
{
  # All unique combinations
  cmb <- rbind(CJ(L= vl8$ID,
                  R= vl8$ID),
               CJ(L= vl12$ID,
                  R= vl12$ID))
  cmb <- unique(cmb)
  cmb <- merge(cmb,
               motifCounts,
               by.x= "L",
               by.y= "ID",
               allow.cartesian= T)
  cmb <- merge(cmb,
               motifCounts,
               by.x= "R",
               by.y= "ID",
               allow.cartesian= T,
               suffixes= c("_L", "_R"))
  # Scaling
  mat <- t(as.matrix(cmb[, !c("L", "R")]))
  mat <- scale(mat)
  # Clustering
  # grid <- somgrid(5, 5, "hexagonal", toroidal= T) #n25
  grid <- somgrid(4, 4, "hexagonal", toroidal= T)
  som <- som(mat,
             grid)
  saveRDS(som,
          "Rdata/uniq_enh_feat/som_codes_enh_motifs_n16.rds")
  # Codes only
  codes <- as.data.table(t(som$codes[[1]]))
  names(codes) <- paste0("motif_SOM_", seq(codes))
  codes <- cbind(cmb[, .(L, R)], codes)
  saveRDS(codes,
          "Rdata/uniq_enh_feat/som_codes_enh_motifs_codesOnly_n16.rds")
}else
  codes <- readRDS("Rdata/uniq_enh_feat/som_codes_enh_motifs_codesOnly_n16.rds")

#---------------------------------------------------------------#
# Gene assignment
#---------------------------------------------------------------#
# Closest promoter
if(!file.exists("Rdata/uniq_enh_feat/closestTss.rds"))
{
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
  saveRDS(closestTss,
          "Rdata/uniq_enh_feat/closestTss.rds")
}else
  closestTss <- readRDS("Rdata/uniq_enh_feat/closestTss.rds")

#---------------------------------------------------------------#
# Chromatin Features
#---------------------------------------------------------------#
if(!file.exists("Rdata/uniq_enh_feat/chromatinEnrich.rds"))
{
  if(!file.exists("db/bed/GSE119708_ATAC_merged.bed"))
  {
    ATAC <- list.files("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/",
                       "ATAC_rep.*.bed$", 
                       full.names = T)
    system(paste(c("cat", ATAC, ">", "db/bed/GSE119708_ATAC_merged.bed"), collapse = " "))
    vl_exportBed(vl_shuffleBed("db/bed/GSE119708_ATAC_merged.bed"), "db/bed/GSE119708_ATAC_shuffled.bed")
  }
  chrom <- data.table(ID= lib$ID)
  chrom$ATAC_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                center = "center", 
                                                                upstream = 500, 
                                                                downstream = 500, 
                                                                ignore.strand = T),
                                         ChIP_bed = "db/bed/GSE119708_ATAC_merged.bed",
                                         Input_bed = "db/bed/GSE119708_ATAC_shuffled.bed")$signalValue)
  
  chrom$H3K27ac_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                   center = "center", 
                                                                   upstream = 500, 
                                                                   downstream = 500, 
                                                                   ignore.strand = T),
                                            ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep2_uniq.bed"),
                                            Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                          "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"))$signalValue)
  
  chrom$H3K4me1_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                   center = "center", 
                                                                   upstream = 500, 
                                                                   downstream = 500, 
                                                                   ignore.strand = T),
                                            ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep2_uniq.bed"),
                                            Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                          "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"))$signalValue)
  
  chrom$H3K4me3_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                   center = "center", 
                                                                   upstream = 500, 
                                                                   downstream = 500, 
                                                                   ignore.strand = T),
                                            ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep2_uniq.bed"),
                                            Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep1_uniq.bed",
                                                          "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep2_uniq.bed"))$signalValue)
  
  chrom$H3K27me3_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                    center = "center", 
                                                                    upstream = 500, 
                                                                    downstream = 500, 
                                                                    ignore.strand = T),
                                             ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27me3_rep1_uniq.bed"),
                                             Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed"))$signalValue)
  
  chrom$Pol2_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                center = "center", 
                                                                upstream = 500, 
                                                                downstream = 500, 
                                                                ignore.strand = T),
                                         ChIP_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_RNAPolII_rep1_uniq.bed",
                                         Input_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed")$signalValue)
  
  chrom$GAF_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                               center = "center", 
                                                               upstream = 500, 
                                                               downstream = 500, 
                                                               ignore.strand = T),
                                        ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep2_uniq.bed"),
                                        Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep1_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep2_uniq.bed"))$signalValue)
  
  chrom$SUHW_log2FC <- log2(vl_enrichBed(regions = vl_resizeBed(lib, 
                                                                center = "center", 
                                                                upstream = 500, 
                                                                downstream = 500, 
                                                                ignore.strand = T),
                                         ChIP_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_SuHw_rep1_uniq.bed",
                                         Input_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_input_rep1_uniq.bed")$signalValue)
  saveRDS(chrom,
          "Rdata/uniq_enh_feat/chromatinEnrich.rds")
}else
  chrom <- readRDS("Rdata/uniq_enh_feat/chromatinEnrich.rds")

#--------------------------------------------------------------------#
# Deep STARR prediction
#--------------------------------------------------------------------#
# Script to run this: 
"/groups/stark/almeida/Projects/Vincent_pairedSTARRseq/results/DeepSTARR_predictions/Run_DeepSTARR.sh"

deep <- fread("/groups/stark/almeida/Projects/Vincent_pairedSTARRseq/results/DeepSTARR_predictions/vl_sequences.fa_predictions_Model_continuous_CNN_pool_120_7_60_3_60_3_pooling_2denses_dropout0.3_valloss2.73.txt")
setnames(deep, c("ID", "seq", "deep_dev", "deep_hk"))
saveRDS(deep[, .(ID, deep_dev, deep_hk)], 
        file = "Rdata/uniq_enh_feat/DeepSTARR_prediction.rds")
#--------------------------------------------------------------------#
# Final object
#--------------------------------------------------------------------#
saveRDS(lib[, .(ID, group, detail, dev_log2FC_STARR200, hk_log2FC_STARR200, dev_log2FC_TWIST, hk_log2FC_TWIST, col)],
        "Rdata/uniq_enh_feat/lib_features.rds")
saveRDS(lib[, .(ID, seqnames, start, end, strand, linker_ID, oligo_full_sequence)],
        "Rdata/uniq_enh_feat/lib_genomic_dat.rds")
saveRDS(data.table(ID= lib$ID, motif_cluster= km_cl),
        "Rdata/uniq_enh_feat/kmCluster_ID.rds")
saveRDS(setnames(data.table(ID= lib$ID, topMotifCounts), colnames(topMotifCounts), function(x) paste0("top_motif_", x)),
        "Rdata/uniq_enh_feat/topMotifKmClusters_IDs.rds")

obj <- list(genomic= "Rdata/uniq_enh_feat/lib_genomic_dat.rds",
            lib= "Rdata/uniq_enh_feat/lib_features.rds",
            motif_cluster= "Rdata/uniq_enh_feat/kmCluster_ID.rds",
            motif_cluster_enrich= "Rdata/uniq_enh_feat/motifEnrichKmClusters.rds",
            top_motifs= "Rdata/uniq_enh_feat/topMotifKmClusters_IDs.rds",
            pairs_motif_codes= "Rdata/uniq_enh_feat/som_codes_enh_motifs_codesOnly_n16.rds",
            closest_TSSs= "Rdata/uniq_enh_feat/closestTss.rds",
            chromatin_features= "Rdata/uniq_enh_feat/chromatinEnrich.rds",
            deepSTARR= "Rdata/uniq_enh_feat/DeepSTARR_prediction.rds",
            add_feature= function(DT, feature) {
              if(!all(c("L", "R") %in% names(DT)) | !is.data.table(DT))
                stop("DT should be a data.table containing L and R columns corresponding to left and right enhancers")
              if(!file.exists(feature))
                stop("feature file not found") else
                  feature <- readRDS(feature)
              if(all(c("L", "R") %in% names(feature)))
                on_cols <- c("L", "R") else if("ID" %in% names(feature))
                  on_cols <- "ID" else
                    stop("feature should contain either L and R columns or a unique ID column corresponding to left/right enhancers")
              # add_cols <- setdiff(names(feature), on_cols)
              
              # Select vars to add to DT
              if(identical(on_cols, c("L", "R")))
              {
                DT <- merge(DT,
                            feature, 
                            by= on_cols,
                            all.x= T)
              }else{
                DT <- merge(DT, 
                            feature, 
                            by.x= "L", 
                            by.y= "ID", 
                            all.x= T)
                DT <- merge(DT, 
                            feature, 
                            by.x= "R", 
                            by.y= "ID", 
                            all.x= T,
                            suffixes= c("_L", "_R"))
              }
                return(DT)
            })

saveRDS(obj, "Rdata/final_300bp_enhancer_features.rds")
