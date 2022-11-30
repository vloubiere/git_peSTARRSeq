setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_dev_pairs_vllib002_with_predictions.rds")

# Select variables of interest
top <- readRDS("Rdata/top_enrich_motifs_residuals_density.rds")
top <- rbind(cbind(top, top), 
             rbind(
               cbind(top[name=="Ebox/CATATG/twi"], top[name=="Trl/1"]), 
               cbind(top[name=="Trl/1"], top[name=="Ebox/CATATG/twi"])))
cols <- c("name1", "name2")
setnames(top, c("name1", "var1", "name2", "var2"))
top[, (cols):= lapply(.SD, as.character), .SDcols= cols]
top[, name:= ifelse(name1==name2, name1, paste0(name1, " x ", name2))]
top$name1 <- top$name2 <- NULL

# Compute motif counts
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
lib <- lib[ID_vl %in% dat[, c(L, R)]]
sel <- unique(top$var1)
counts <- vl_motif_counts(lib$enh_sequence, 
                          sel= sel, 
                          genome= "dm3", 
                          collapse_overlapping = T)
counts <- cbind(lib, counts)

# Plotting object
pl <- top[, (dat[, .(L, R, residuals, indL, indR)]), (top)]
pl[, countsL:= counts[L, var1, on= "ID_vl", with= F][[1]], var1]
pl[, countsR:= counts[R, var2, on= "ID_vl", with= F][[1]], var2]
# Counts cutoff
pl[, minL:= min(c(4, quantile(countsL, 0.99))), name]
pl[, minR:= min(c(4, quantile(countsR, 0.99))), name]
pl[countsL>minL, countsL:= minL]
pl[countsR>minR, countsR:= minR]
pl[, cut:= paste0(countsL, "_", countsR)]
pl[, N:= .N, .(name, cut)]
pl <- pl[(N>=50)]

#--------------------------------#
# PLOT
#--------------------------------#
pdf("pdf/draft/heatmap_motif_residuals_prediction.pdf", 4, 4)
par(mar= c(4, 4, 4, 4),
    tcl= -0.2,
    las= 1)
pl[, {
  .lm <- summary(lm(residuals~indL*indR+cut))$coefficients
  .lm <- as.data.table(.lm, keep.rownames = T)[grepl("^cut", rn)]
  .lm[, c("L", "R") := lapply(tstrsplit(gsub("^cut", "", rn), "_"), as.integer)]

  mat <- dcast(.lm, L~R, value.var = "Estimate")
  mat <- as.matrix(mat, 1)
  mat <- mat[nrow(mat):1,]
  
  vl_heatmap(mat,
             cluster_rows = F,
             cluster_cols = F,
             breaks= c(-1,0,1),
             main= name,
             show_legend= F)
  vl_heatkey(breaks = c(-1,0,1),
             left= par("usr")[2]+strwidth("M"),
             col= c("cornflowerblue", "white", "red"),
             main = "lm coeff.")
  print(as.character(name))
  .SD
}, name]
dev.off()
