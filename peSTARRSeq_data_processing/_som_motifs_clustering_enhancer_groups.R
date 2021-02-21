sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

addClustersToLib <- function(lib)
{
  dat <- copy(lib)
  
  #-----------------------#
  # Chromatin clustering
  #-----------------------#
  chrom <- dat[, .(ID, ATAC_log2FC, H3K4me1_log2FC, H3K4me3_log2FC, H3K27ac_log2FC, H3K27me3_log2FC)]
  chrom <- na.omit(chrom)
  # Clip outliers
  cols <- grep("log2FC$", colnames(chrom), value = T)
  chrom[, (cols) := lapply(.SD, function(x) {
    lims <- quantile(x, c(0.05, 0.95))
    x[x<lims[1]] <- lims[1]
    x[x>lims[2]] <- lims[2]
    return(x)
  }), .SDcols= cols]
  chrom <- melt.data.table(chrom, id.vars= "ID")
  
  pdf("pdf/peSTARRSeq/chrom_enhancers_clustering.pdf")
  chrom_cl <- my_heatmap(chrom, row.BY = "ID", col.BY = "variable", value.var = "value", cutree_rows = 7)
  dev.off()
  
  saveRDS(chrom_cl, "Rdata/processed_peSTARRSeq_data/chromatin_enhancers_clustering.rds")
  
  #-----------------------#
  # Motifs clustering
  #-----------------------#
  motifs <- dat[, motif__cisbp__M2330:motif__predrem__nrMotif1321]
  motifs[, ID:= dat$ID]
  motifs <- na.omit(motifs[!grepl("NegativeRegions|^control", ID)])
  
  # Clip outliers
  cols <- grep("^motif", colnames(motifs), value = T)
  motifs[, (cols) := lapply(.SD, scale), .SDcols= cols]
  motifs <- melt.data.table(motifs, id.vars= "ID")
  
  pdf("pdf/peSTARRSeq/motifs_enhancers_clustering.pdf", width = 18, height = 12)
  layout(matrix(1:2, nrow= 1), widths = c(1, 0.3))
  par(mar= c(21,17.5,5,7))
  mot_cl <- my_heatmap(motifs, row.BY = "ID", col.BY = "variable", value.var = "value", 
                       cutree_rows = 20, clustering_dist_rows = "pearson", breaks = c(-3, 0, 5), plot_dendro_row = F)
  mot_cl[dat, group:= i.group, on= "ID"]
  groups <- dcast(mot_cl, ycoor~group, fun.aggregate = function(x) ifelse(any(x>0), 1, 0))
  groups <- as.matrix(groups[, -1], 1)
  groups <- groups[nrow(groups):1,]
  par(mar= c(21,2,5,7))
  my_pheatmap(groups, cluster_rows = F, cluster_cols= F)
  dev.off()
  
  saveRDS(mot_cl, "Rdata/processed_peSTARRSeq_data/motifs_enhancers_clustering.rds")
  
  # Add to lib
  lib[chrom_cl, chromatin_cl:= i.row_cl, on= "ID"]
  lib[mot_cl, motif_cl:= i.row_cl, on= "ID"]
  
  return(lib)
}
