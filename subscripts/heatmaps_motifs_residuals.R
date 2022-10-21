setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- unique(dat[actClassL!= "inactive", .(L, indL)])[indL<=quantile(indL, 0.9), L]
QR <- unique(dat[actClassR!= "inactive", .(R, indR)])[indR<=quantile(indR, 0.9), R]
dat <- dat[L %in% QL & R %in% QR]

# Select variables of interest
top <- readRDS("Rdata/top_enrich_motifs_residuals_density.rds")
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
lib <- lib[ID_vl %in% dat[, c(L, R)]]
counts <- vl_motif_counts(lib$enh_sequence, sel= top$variable)
counts <- cbind(lib, counts)

#--------------------------------#
# PLOT
#--------------------------------#
pdf("pdf/draft/heatmap_motif_residuals_prediction.pdf", 4, 4)
par(mar= c(4, 4, 4, 4),
    tcl= -0.2,
    las= 1)
top[, {
  sub <- data.table(dat[, .(residuals, indL, indR)],
                    L= counts[dat, variable, on= "ID_vl==L", with= F][[1]],
                    R= counts[dat, variable, on= "ID_vl==R", with= F][[1]])
  QL <- floor(quantile(sub$L, 0.99))
  QR <- floor(quantile(sub$R, 0.99))
  sub[L>QL, L:= QL]
  sub[R>QR, R:= QR]
  sub[, cut:= paste0(L, "_", R)]
  .lm <- summary(lm(residuals~indL*indR+cut, sub))$coefficients
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
}, .(name, variable)]
dev.off()
