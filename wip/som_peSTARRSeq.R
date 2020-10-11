setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[!is.na(median_L) & !is.na(median_R) & !is.na(log2FC_add) & !is.na(log2FoldChange) & !is.na(diff) & enh_L!=enh_R]
feat <- readRDS("Rdata/library/lib_features.rds")

#------------------------------------------------------------#
# 1- Train SOM
#------------------------------------------------------------#

if(!file.exists("Rdata/motifs/som_peSTARR_supervised.rds"))
{
  # obj <- list(act= as.matrix(dat[, .(ID, median_L, median_R, log2FC_add, log2FoldChange, diff)], 1),
  #             chrom_L= as.matrix(merge(dat[, .(ID= enh_L)], feat[, .SD, .SDcols= patterns("ID|log2FC$")]), 1),
  #             chrom_R= as.matrix(merge(dat[, .(ID= enh_R)], feat[, .SD, .SDcols= patterns("ID|log2FC$")]), 1),
  #             motif_L= as.matrix(merge(dat[, .(ID= enh_L)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1),
  #             motif_R= as.matrix(merge(dat[, .(ID= enh_R)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1),
  #             group_L= factor(merge(dat[, .(ID= enh_L)], feat[, .(ID, group)])$group),
  #             group_R= factor(merge(dat[, .(ID= enh_R)], feat[, .(ID, group)])$group))
  
  obj <- list(coop= as.matrix(dat[, .(ID, diff)], 1),
              act= as.matrix(dat[, .(ID, median_L, median_R, log2FoldChange)], 1),
              motif_L= as.matrix(merge(dat[, .(ID= enh_L)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1),
              motif_R= as.matrix(merge(dat[, .(ID= enh_R)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1))
  
  # Supervised
  mygrid <- somgrid(xdim= 18, ydim= 18, topo = 'hexagonal', toroidal = T)
  set.seed(1234)
  # som.model <- supersom(obj, user.weights = c(1,3,3), grid = mygrid, rlen = 100, maxNA.fraction = 0.3)
  som.model <- supersom(obj, grid = mygrid, rlen = 100, maxNA.fraction = 0.3)
  saveRDS(som.model, "Rdata/motifs/som_peSTARR_supervised_3.rds")
}

#------------------------------------------------------------#
# 2- plot SOMs
#------------------------------------------------------------#
# Import codes
som <- readRDS("Rdata/motifs/som_peSTARR_supervised_3.rds")
dmat <- as.data.table(som$codes, keep.rownames = F)
codes <- melt(dmat, measure.vars = colnames(dmat))
# Kmeans clustering
kcl <- kmeans(dmat[, coop.diff], 6)$cluster

pdf("pdf/peSTARRSeq/som_peSTARRFC_vs_motifs.pdf", 60, 60)
par(mfrow=c(8,9))
plot(som, "property", property= kcl, palette.name= rainbow, shape= "straight", border= NA, main= "clusters")
Cc <- colorRampPalette(c("black", "blue", "yellow", "tomato"))
codes[, 
      {
        plot(som, "property", property= value, palette.name= Cc, shape= "straight", border= NA, main= variable)
        add.cluster.boundaries(som, kcl, col="white")
      }, variable]
dev.off()
