setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(GenomicRanges)
require(motifmatchr)
require(seqLogo)
require(TFBSTools)
require(kohonen)

lib <- readRDS("Rdata/uniq_library_final.rds")
lib <- lib[!detail=="ecoli"]

#-----------------------------------------------#
# 1- Motif counts
#-----------------------------------------------#
# Count hits low cutoff
if(!file.exists("Rdata/all_motifs_counts_lib.rds"))
{
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  DT <- as.data.table(TF_clusters_PWMs$metadata)
  hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds, GRanges(lib$coor), genome= "dm3", p.cutoff= 5e-4, bg="even", out= "scores")
  counts <- as.matrix(motifCounts(hit))
  colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds)
  rownames(counts) <- lib$ID
  counts <- as.data.table(counts, keep.rownames = T)
  saveRDS(counts, "Rdata/all_motifs_counts_lib.rds")
}

#-----------------------------------------------#
# 2- Select motifs enriched in at least one of the lib groups
#-----------------------------------------------#
if(!file.exists("Rdata/enriched_motifs_counts.rds"))
{
  counts <- readRDS("Rdata/all_motifs_counts_lib.rds")
  
  dat <- melt(counts, id.vars= "rn")
  dat[lib, group:= i.group, on= "rn==ID"]
  sel <- dat[group!="control", any(value>0), .(variable, group)][(V1), .(variable, group)]
  setkey(sel, group)
  
  mot <- list()
  for(c_gr in c("dev", "hk", "inducible", "OSC"))
  {
    sub <- dat[group %in% c(c_gr, "control") & variable %in% sel[c_gr, variable]]
    sub <- sub[, fisher.test(table(group!="control", value>0))[c("estimate", "p.value")], variable]
    mot[[c_gr]] <- sub[estimate>2 & p.value<0.00001, variable]
  }
  enriched_motifs_counts <- dat[variable %in% unique(unlist(mot))]
  saveRDS(enriched_motifs_counts, "Rdata/enriched_motifs_counts.rds")
}

#-----------------------------------------------#
# 3- SOM motifs
#-----------------------------------------------#
if(!file.exists("Rdata/som_enriched_motifs.rds"))
{
  dat <- readRDS("Rdata/enriched_motifs_counts.rds")
  dmat <- dcast(dat, variable~rn, value.var = "value")
  mat <- log2(as.matrix(dmat, 1)+1)
  mat <- scale(mat) 
  mat <- apply(mat, 2, function(x)
  {
    lim <- quantile(x, 0.025, 0.975)
    x[x<lim[1]] <- lim[1]
    x[x>lim[2]] <- lim[2]
    return(x)
  })
  
  mygrid <- somgrid(xdim= 6, ydim= 6, topo = 'hexagonal', toroidal = T)
  set.seed(1234)
  som.model <- supersom(mat, grid = mygrid, rlen = 500)
  saveRDS(som.model, "Rdata/som_enriched_motifs.rds")
  
  pdf("pdf/som_enriched_motifs_clustering_diag_plots.pdf")
  par(mfrow= c(4, 4))
  for(type in c("changes", "dist.neighbours", "counts", "mapping", "quality"))
  {
    plot(som.model, type= type, palette.name = colorRampPalette(c("blue", "yellow")), shape= "straight", border= "black")
  }
  dev.off()
  
  pdf("pdf/som_enriched_motifs_clustering_pheatmap.pdf")
  pl <- mat[order(som.model$unit.classif),]
  my_pheatmap(pl, cluster_rows = F, col= colorRampPalette(c("blue", "yellow"))(100), plot_dendro_col = F)
  abline(h= 1-which(diff(som.model$unit.classif[order(som.model$unit.classif)])!=0)/nrow(mat), lwd= 1, col= "white")
  dev.off()
}


















mat <- mat[order(som.model$unit.classif),]
pdf("test.pdf")
my_pheatmap(mat, cluster_rows = F, col= colorRampPalette(c("blue", "yellow"))(100))
abline(h= 1-which(diff(som.model$unit.classif[order(som.model$unit.classif)])!=0)/791, lwd=2, col= "white")
dev.off()

## use hierarchical clustering to cluster the codebook vectors
fviz_nbclust(as.matrix(object.distances(som.model, "codes")), FUN = hcut, method = "wss", k.max = 30)
som.model$hclust <- cutree(hclust(object.distances(som.model, "codes")), 9)


# diagnostic plots
Cc <- colorRampPalette(c("black", "blue", "yellow"))

par(mfrow= c(4, 4))
for(type in c("changes", "dist.neighbours", "counts", "mapping", "quality"))
{
  plot(som.model, type= type, palette.name = Cc, shape= "straight", border= "black")
  # if(type != "changes")
  # {
  #   add.cluster.boundaries(som.model, som.model$hclust, lwd= 3, col= "red")
  # }
}
for(what in names(som.model$codes))
{
  plot(som.model, type= "codes", whatmap= what, shape= "straight", border= "black")
  add.cluster.boundaries(som.model, som.model$hclust, lwd= 3, col= "red")
}
