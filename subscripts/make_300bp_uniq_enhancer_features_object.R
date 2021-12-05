setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(motifmatchr)

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
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
sel <- TF_clusters_PWMs$metadata$motif_name[!is.na(TF_clusters_PWMs$metadata$FBgn)]
sel <- name(TF_clusters_PWMs$All_pwms_log_odds) %in% sel
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], 
                   DNAStringSet(lib$oligo_full_sequence),
                   p.cutoff= 5e-4, 
                   bg="even", 
                   out= "scores")
counts <- as.matrix(motifCounts(hit))
colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds)[sel]
rownames(counts) <- lib$ID
# Select motifs enriched in at least one of the lib groups
.m <- melt(as.data.table(counts, keep.rownames = T), 
           id.vars = "rn")
.m[lib, "group":= i.group, on= "rn==ID"]
.m[as.data.table(TF_clusters_PWMs$metadata[, c("motif_name", "FBgn")]), "FBgn":= i.FBgn, on= "variable==motif_name"]
cls <- data.table(cl= unique(.m[group!="control", group]))
.m <- .m[, cls[, {
  .tab <- table(group==cl, value>1)
  if(identical(dim(.tab), c(2L,2L)))
  {
    .c<- fisher.test(.tab)
    .(OR= .c$estimate, pval= .c$p.value)
  }else
    .(OR= as.numeric(NA), pval= as.numeric(NA))
}, cl], .(variable, FBgn)]
.m[, padj:= p.adjust(pval, "fdr")]
enriched_motifs <- unique(.m[padj<0.001 & OR>1, variable])
saveRDS(.m, "Rdata/motif_enrichment_per_group.rds")
# Cluster enhancers based on enriched motifs
enr <- counts[, colnames(counts) %in% enriched_motifs]
enr <- log2(enr+1)
# Clip outliers
enr <- apply(enr, 2, function(x){
  lim <- quantile(x, c(0.05, 0.95))
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  return(x)
})
# Cluster and add to lib objectsave heatmap clustering
final_cl <- vl_heatmap(enr, 
                       cutree_rows = 12,
                       clustering_distance_cols = "spearman",
                       plot = F)
lib[final_cl$result, motif_cl:= i.rn_cl, on= "ID==rn"]
# Plot clustering diag
pdf("pdf/motifs_clustering/heatmap_motif_clustering_3200bp_uniq_enhancers.pdf", 
    height = 10)
plot(final_cl,
     breaks = c(0, 3), 
     col = c("blue", "yellow"), 
     cutree_rows= 12,
     cutree_cols = 12)
association <- CJ(group= unique(lib$group),
                  motif_cl= unique(lib$motif_cl))
association[, c("OR", "pval"):= 
              fisher.test(table(lib$group==group, lib$motif_cl==motif_cl))[c("estimate", "p.value")], .(group, motif_cl)]
mat <- dcast(association, group~motif_cl, value.var = "OR")
pval <- dcast(association, group~motif_cl, value.var = "pval")
pval <- as.matrix(pval, 1)
mat <- as.matrix(mat, 1)
mat <- log2(mat+min(mat[mat>0]))
mat[pval>0.05] <- NA
vl_heatmap(mat, 
           breaks = c(-2, 0, 2), 
           cluster_rows=F, 
           cluster_cols=F)
dev.off()

#---------------------------------------------------------------#
# Chromatin features
#---------------------------------------------------------------#
# ChIP-Seq enrichment
ChIP <- lib[, .(ID= .(ID)), .(seqnames, start, end)]
ChIP$ATAC_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep1_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep2_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep3_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep4_uniq.bed"),
                                         peaks = ChIP[, .(seqnames, start, end)], 
                                         ext_peaks = 1000)$log2_enr
ChIP$H3K27ac_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep2_uniq.bed"),
                                            Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                          "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"),
                                            peaks = ChIP[, .(seqnames, start, end)], 
                                            ext_peaks = 1000)$log2_enr
ChIP$H3K4me1_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep2_uniq.bed"),
                                           Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"),
                                           peaks = ChIP[, .(seqnames, start, end)], 
                                           ext_peaks = 1000)$log2_enr
ChIP$H3K4me3_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep1_uniq.bed",
                                                         "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep2_uniq.bed"),
                                            Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep1_uniq.bed",
                                                          "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep2_uniq.bed"),
                                            peaks = ChIP[, .(seqnames, start, end)], 
                                            ext_peaks = 1000)$log2_enr
ChIP$H3K27me3_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27me3_rep1_uniq.bed"),
                                             Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed"),
                                             peaks = ChIP[, .(seqnames, start, end)], 
                                             ext_peaks = 1000)$log2_enr
ChIP$Pol2_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_RNAPolII_rep1_uniq.bed"),
                                         Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed"),
                                         peaks = ChIP[, .(seqnames, start, end)], 
                                         ext_peaks = 1000)$log2_enr
ChIP$GAF_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep1_uniq.bed",
                                                     "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_GAF_rep2_uniq.bed"),
                                        Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep1_uniq.bed",
                                                      "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE40646_input_rep2_uniq.bed"),
                                        peaks = ChIP[, .(seqnames, start, end)], 
                                        ext_peaks = 1000)$log2_enr
ChIP$SUHW_log2FC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_SuHw_rep1_uniq.bed"),
                                         peaks = ChIP[, .(seqnames, start, end)], 
                                         ext_peaks = 1000)$log2_enr
ChIP <- ChIP[, .(ID= unlist(ID)), setdiff(names(ChIP), "ID")]
lib <- lib[ChIP[, ATAC_log2FC:ID], on= "ID"]

#------------------------------------------------------------------------------------------------#
# SAVE
#------------------------------------------------------------------------------------------------#
saveRDS(lib, "Rdata/final_300bp_enhancer_features.rds")
