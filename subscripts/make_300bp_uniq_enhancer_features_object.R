setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(motifmatchr)
require(rtracklayer)
require(kohonen)

#------------------------------------------------------------------------------------------------------------------------------#
# Make enhancers unique and recompile features that were used for design 
#------------------------------------------------------------------------------------------------------------------------------#
# Merge BA and my TWIST into a single non-redundant clean object
vl8 <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist008_112019.rds"))
names(vl8)[names(vl8)=="ID_vl"] <- "ID"
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
# Collapse the two libraries (some oligos are shared)
lib <- unique(lib)
# Add values from BA screen
BA <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/BA_300bp_TWIST_STARRSeq.txt")
lib[, c("dev_log2FoldChange", "hk_log2FoldChange"):= 
      BA[.BY, .(dev_log2FoldChange, hk_log2FoldChange), 
         on= c("seqnames", "start", "end", "strand")], .(seqnames, start, end, strand)]
# Uniq groups
lib[group=="ecdysone", c("group", "detail"):= .("inducible", "ecdysone")]
lib[group=="heatshock", c("group", "detail"):= .("inducible", "heatshock")]
lib[group=="Repressor", group:= "repressor"]
lib[BA[enhancer_group=="shared"], group:= "shared", on= c("seqnames", "start", "end", "strand")]
lib[, group:= factor(group, 
                     levels= c("repressor", "SUHW_peak", "inducible", "OSC", "CP", "control", "DHS_peak", "hk", "shared", "dev"))]
# Uniq details
lib[group== "shared", 
    detail:= cut(rowMeans(.SD), 
                 quantile(rowMeans(.SD), seq(0, 1, length.out= 5)), 
                 labels = c("inactive", "weak", "medium", "strong"), 
                 include.lowest= T),
    .SDcols= c("hk_log2FoldChange", "dev_log2FoldChange")]
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
class_Cc <- data.table(group= c("hk", "dev", "shared", "OSC", "inducible", "control", "CP", "DHS_peak", "repressor", "SUHW_peak"),
                       col= c("tomato", "#74C27A", "darkgreen", "black", "gold", "lightgrey", "cyan3", "royalblue2", "deeppink2", "darkorchid3"))
lib <- lib[class_Cc, on= "group"]

#------------------------------------------------------------------------------------------------------------------------------#
# TF motifs
#------------------------------------------------------------------------------------------------------------------------------#
# Count hits low cutoff
hits <- vl_motif_counts(sequences =  lib$oligo_full_sequence)
# Select motifs enriched in at least one of the lib groups
enr <- vl_motif_cl_enrich(hits, 
                          cl_IDs = lib$group, 
                          plot = F, 
                          N_top = 10)
motifs <- as.character(unique(enr[padj<0.001, motif]))
motifCounts <- hits[, ..motifs]
motifCounts[, ID:= lib$ID]
setcolorder(motifCounts, "ID")
# SOM clustering (~100Go of RAM!!)
if(!file.exists("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs_codesOnly_n16.rds"))
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
          "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs_n16.rds")
  # Codes only
  codes <- as.data.table(t(som$codes[[1]]))
  names(codes) <- paste0("motif_SOM_", seq(codes))
  codes <- cbind(cmb[, .(L, R)], codes)
  saveRDS(codes,
          "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs_codesOnly_n16.rds")
}else
  codes <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs_codesOnly_n16.rds")

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

#---------------------------------------------------------------#
# Chromatin Features
#---------------------------------------------------------------#
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

#--------------------------------------------------------------------#
# Final object
#--------------------------------------------------------------------#
obj <- list(genomic= lib[, .(ID, seqnames, start, end, strand, linker_ID, oligo_full_sequence)],
            lib= lib[, .(ID, group, detail, dev_log2FoldChange, hk_log2FoldChange, col)],
            top_motifs= cbind(ID= lib$ID, hits[,  ..motifs]),
            motif_padj= enr[motif %in% motifs],
            pairs_motif_codes= codes,
            closest_TSSs= closestTss,
            chromatin_features= chrom,
            add_feature= function(DT, feature) {
              if(!all(c("L", "R") %in% names(DT)) | !is.data.table(DT))
                stop("DT should be a data.table containing L and R columns corresponding to left and right enhancers")
              
              # Select vars to add to DT
              if(all(c("L", "R") %in% names(feature)))
              {
                return(DT[feature, on= c("L", "R"), nomatch= NULL])
              }else
              {
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
                return(DT)
              }
            })
saveRDS(obj, "Rdata/final_300bp_enhancer_features.rds")
