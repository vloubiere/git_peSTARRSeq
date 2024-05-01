setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
S2 <- readRDS("db/linear_models/lasso_DSCP_large_WT_residuals.rds")
ECD <- readRDS("db/linear_models/lasso_DSCP_ECD_WT_residuals.rds")
OSC <- readRDS("db/linear_models/lasso_DSCP_OSC_WT_residuals.rds")
dat <- list("S2 - ecd."= S2,
            "S2 + ecd."= ECD,
            "OSC cells"= OSC)
dat <- lapply(dat, function(x) as.data.table(as.matrix(x$log2FoldChange$beta), keep.rownames= T))
dat <- rbindlist(dat, idcol = "cdition")
dat[, cdition:= factor(cdition, 
                       c("S2 - ecd.",
                         "S2 + ecd.",
                         "OSC cells"))]

# Select top motifs ----
dat[, rank:= rank(-s0), .(cdition, grepl("__L", rn))]
sel <- dat[rank<=3, rn]
sel <- gsub("__L$|__R$", "", sel)
sel <- unique(sel)
sel <- paste0(sel, collapse= "|")
dat <- dat[grepl(sel, rn)]

# Dcast Left and right motif ----
mat <- dcast(dat, rn~cdition, value.var = "s0")
mat <- as.matrix(mat, 1)
matL <- mat[grepl("__L$", rownames(mat)),]
matR <- mat[grepl("__R$", rownames(mat)),]

# order based on hclust ----
ord <- hclust(dist(matR))$order
matL <- matL[ord,]
matR <- matR[ord,]

# Retrieve logos ----
IDs <- readRDS("db/motif_counts/lib8_motifs_IDs.rds")
pwms <- IDs[gsub("^.*?__|__[^__]*$", "", rownames(matL)), pwm, on= "motif_ID"]
rownames(matL) <- unlist(tstrsplit(rownames(matL), "__", keep=1))
rownames(matR) <- unlist(tstrsplit(rownames(matR), "__", keep=1))

pdf("pdf/draft/LASSO_coeff_per_cell_type.pdf", 3.9, 3.8)
vl_par(mai= c(.9,.05,.9,0),
       mfrow= c(1,2),
       omi= c(0, 1.8, 0, 0.9),
       font.main= 1,
       cex.main= 7/12,
       xpd= NA,
       bty= "o")
vl_heatmap(matL,
           cluster.rows = F,
           cluster.cols = F,
           breaks = c(-.4, 0, .4),
           show.legend = F,
           main= "5' enhancer",
           tilt.colnames = T,
           box.lwd = .5)
vl_seqlogo(pwms,
           x= par("usr")[1]-strwidth(rownames(matL), cex= par("cex.lab"))-strwidth("M"),
           y= rev(seq(nrow(matL))),
           pos = 2,
           add= T,
           cex.width = .5,
           min_content = .1)
vl_heatmap(matR,
           cluster.rows = F,
           cluster.cols = F,
           breaks = c(-.4, 0, .4),
           legend.title = "LASSO coef.",
           show.rownames = F,
           main= "3' enhancer",
           tilt.colnames = T, 
           legend.cex = 6/12,
           box.lwd = .5)
dev.off()

# comb <- merge(dat$S2, dat$ECD, by= "rn", suffixes= c("_S2", "_ECD"))
# comb <- merge(comb, dat$OSC, by= "rn")
# 
# 
# 
# setnames(comb, "s0", "s0_OSC")
# mat <- as.matrix(comb, 1)
# pca <- prcomp(mat)
# tab <- pca$x
# 
# plot(tab[,"PC1"], tab[,"PC2"])
# text(tab[,"PC1"], tab[,"PC2"], unlist(tstrsplit(rownames(tab), "__", keep= 1)))
# 
# dat <- merge(as.data.table(as.matrix(mod$residuals$beta), keep.rownames= T),
#              as.data.table(as.matrix(mod$log2FoldChange$beta), keep.rownames= T),
#              by= "rn",
#              suffixes= c("res", "act"))
# 
# plot(dat$s0res, dat$s0act)
# text(dat$s0res, dat$s0act, unlist(tstrsplit(dat$rn, "__", keep= 1)))
# 
# # test <- readRDS("db/linear_models/FC_lasso_DSCP_ECD_WT_residuals_predictions.rds")
# # test <- readRDS("db/linear_models/FC_lasso_DSCP_large_WT_residuals_predictions.rds")
# test <- readRDS("db/linear_models/FC_lasso_DSCP_OSC_WT_residuals_predictions.rds")
# smoothScatter(test$log2FoldChange, test$predicted_lasso)
# vl_model_eval(test$log2FoldChange, test$predicted_lasso)
# smoothScatter(test$residuals, test$predicted_residuals_lasso)
# vl_model_eval(test$residuals, test$predicted_residuals_lasso)
