setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[!is.na(median_L) & !is.na(median_R) & !is.na(log2FC_add) & !is.na(log2FoldChange) & !is.na(diff) & enh_L!=enh_R]
feat <- readRDS("Rdata/library/lib_features.rds")

# obj <- list(act= as.matrix(dat[, .(ID, median_L, median_R, log2FC_add, log2FoldChange, diff)], 1),
#             chrom_L= as.matrix(merge(dat[, .(ID= enh_L)], feat[, .SD, .SDcols= patterns("ID|log2FC$")]), 1),
#             chrom_R= as.matrix(merge(dat[, .(ID= enh_R)], feat[, .SD, .SDcols= patterns("ID|log2FC$")]), 1),
#             motif_L= as.matrix(merge(dat[, .(ID= enh_L)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1),
#             motif_R= as.matrix(merge(dat[, .(ID= enh_R)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1),
#             group_L= factor(merge(dat[, .(ID= enh_L)], feat[, .(ID, group)])$group),
#             group_R= factor(merge(dat[, .(ID= enh_R)], feat[, .(ID, group)])$group))

obj <- list(act= as.matrix(dat[, .(ID, median_L, median_R, log2FC_add, log2FoldChange, diff)], 1),
            motif_L= as.matrix(merge(dat[, .(ID= enh_L)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1),
            motif_R= as.matrix(merge(dat[, .(ID= enh_R)], feat[, .SD, .SDcols= patterns("ID|^motif")]), 1))

mygrid <- somgrid(xdim= 20, ydim= 20, topo = 'hexagonal', toroidal = T)
if(!file.exists("Rdata/som_peSTARR_supervised.rds"))
{
  # Supervised
  set.seed(1234)
  som.model <- supersom(obj, user.weights = c(1,3,3), grid = mygrid, rlen = 100, maxNA.fraction = 0.3)
  saveRDS(som.model, "Rdata/som_peSTARR_supervised.rds")
}
# if(!file.exists("Rdata/som_peSTARR_unsupervised.rds"))
# {
#   # Unsupervised
#   set.seed(1234)
#   som.model <- supersom(obj, grid = mygrid, rlen = 200, maxNA.fraction = 0.3)
#   saveRDS(som.model, "Rdata/som_peSTARR_unsupervised.rds")
# }
