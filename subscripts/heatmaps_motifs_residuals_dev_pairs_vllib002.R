setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")
dat <- dat[!grepl("control", L) & !grepl("control", R) 
           & (actL!="Inactive" | actR!="Inactive")]

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
pl[, c("cutL", "cutR"):= lapply(.SD, function(x) 
{
  ct <- c(-Inf, 0, unique(quantile(x[x>0], c(0, .25, 0.5, .75))), Inf)
  lab <- ct[-1]
  lab[length(lab)] <- max(x, na.rm= T)
  lab[c(0,diff(lab))>1] <- paste0(lab[diff(lab)>1]+1, "-", lab[c(0,diff(lab))>1])
  .c <- cut(x, ct, lab)
}), name, .SDcols= c("countsL", "countsR")]
pl[, cut:= paste0(cutL, "_", cutR)]

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
  .lm[, c("L", "R") := tstrsplit(gsub("^cut", "", rn), "_")]

  mat <- dcast(.lm, L~R, value.var = "Estimate")
  mat <- as.matrix(mat, 1)
  mat <- mat[nrow(mat):1,]
  pval <- dcast(.lm, L~R, value.var = "Pr(>|t|)")
  pval <- as.matrix(pval, 1)
  pval <- pval[nrow(pval):1,]
  pval <- matrix(cut(pval, 
                     c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf),
                     c("****", "***", "**", "*", "N.S")),
                 ncol= ncol(pval),
                 nrow= nrow(pval))
  pairs <- dcast(.SD, cutL~cutR, fun.aggregate = length)
  pairs <- as.matrix(pairs, 1)
  pairs <- pairs[nrow(pairs):1,]
  pval <- abind::abind(pval, pairs, along = 3)
  pval <- apply(pval, 1:2, function(x) paste0(x, collapse= "\n"))
  
  vl_heatmap(mat,
             cluster_rows = F,
             cluster_cols = F,
             breaks= c(-1,0,1),
             main= name,
             show_legend= F, 
             display_numbers= T,
             display_numbers_matrix= pval, 
             display_numbers_cex= 0.8)
  vl_heatkey(breaks = c(-1,0,1),
             left= par("usr")[2]+strwidth("M"),
             col= c("cornflowerblue", "white", "red"),
             main = "lm coeff.")
  print(as.character(name))
  .SD
}, name]
dev.off()

file.show("pdf/draft/heatmap_motif_residuals_prediction.pdf")
