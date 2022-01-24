require(vlfunctions)
require(motifmatchr)
require(rtracklayer)
require(kohonen)

#---------------------------------------------------------------#
# Recompiles features that were used fot design
#---------------------------------------------------------------#
uniq_defining_features <- function()
{
  #---------------------------------------------------------------#
  # Clean object containing available data
  #---------------------------------------------------------------#
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
  lib[, group:= factor(group, 
                       c("repressor", "SUHW_peak", "inducible", "OSC", "CP", "control", "DHS_peak", "dev", "hk", "shared"))]
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
  class_Cc <- data.table(group= c("hk", "dev", "shared", "OSC", "inducible", "control", "CP", "DHS_peak", "repressor", "SUHW_peak"),
                         col= c("tomato", "#74C27A", "darkgreen", "black", "gold", "lightgrey", "cyan3", "royalblue2", "deeppink2", "darkorchid3"))
  lib <- lib[class_Cc, on= "group"]
}

uniq_add_motif_somCodes <- function(lib, rerun= F)
{
  if(rerun)
  {
    #---------------------------------------------------------------#
    # Seclect motifs used for clusteirng
    #---------------------------------------------------------------#
    # Count hits low cutoff
    hits <- vl_motif_counts(sequences =  lib$oligo_full_sequence)
    # Select motifs enriched in at least one of the lib groups
    enr <- vl_motif_cl_enrich(hits, 
                              cl_IDs = lib$group, 
                              plot = F, 
                              N_top = 10)
    motifs <- as.character(unique(enr[padj<0.001, motif]))
    data <- hits[, ..motifs]
    data[, ID:= lib$ID]
    #---------------------------------------------------------------#
    # Compute enhancer pairs motif combinations
    #---------------------------------------------------------------#
    # Import TWIST libraries
    vl8 <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_twist008_112019.rds"))
    names(vl8)[names(vl8)=="ID_vl"] <- "ID"
    vl12 <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))
    # All combinations unique
    cmb <- rbind(CJ(L= vl8$ID,
                    R= vl8$ID),
                 CJ(L= vl12$ID,
                    R= vl12$ID))
    cmb <- unique(cmb) 
    cmb <- merge(cmb, 
                 data, 
                 by.x= "L",
                 by.y= "ID", 
                 allow.cartesian= T)
    cmb <- merge(cmb, 
                 data, 
                 by.x= "R",
                 by.y= "ID", 
                 allow.cartesian= T, 
                 suffixes= c("_L", "_R"))
    #---------------------------------------------------------------#
    # SOM clustering
    #---------------------------------------------------------------#
    # Scaling
    mat <- t(as.matrix(cmb[, !c("L", "R")]))
    mat <- scale(mat)
    # Cluster and add to lib objectsave heatmap clustering
    grid <- somgrid(5, 5, "hexagonal", toroidal= T)
    som <- som(mat, 
               grid)
    saveRDS(som, 
            "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs.rds")
    codes <- as.data.table(t(som$codes[[1]]))
    names(codes) <- paste0("motif_SOM_", seq(codes))
    codes <- cbind(cmb[, .(L, R)], codes)
    saveRDS(codes,
            "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs_codesOnly.rds")
  }else
  {
    return(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/som_codes_enh_motifs_codesOnly.rds"))
  }
}
