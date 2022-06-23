setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import dataset
if(!exists("feat"))
  feat <- readRDS("Rdata/final_300bp_enhancer_features_w_motifs.rds")
if(!exists("vl_screen"))
  vl_screen <- readRDS("Rdata/final_results_table.rds")

#-----------------------------------------------#
# Train activity based models and compute residuals
#-----------------------------------------------#
clean <- vl_screen[vllib=="vllib002" & class== "enh./enh."]
clean <- clean[L %in% feat[group %in% c("hk", "dev", "shared"), ID] 
               & R %in% feat[group %in% c("hk", "dev", "shared"), ID]]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
clean[, diff:= log2FoldChange-predict(model, new= clean)]

#-----------------------------------------------#
# Format data
#-----------------------------------------------#
# Counts and classes
sub <- feat[ID %in% unique(unlist(clean[, .(L, R)]))]
classes <- sub[, .(ID, group, variable= "class", col)]
counts <- melt(sub[, .(ID, group, Dref_motif, kay_Jra_1_motif, grn_4_motif, grn_24_motif)], id.vars= c("ID", "group"))
# counts[, value:= log2(value+1)]
ind_L <- melt(unique(clean[, .(ID= L, ind_L= median_L)]), id.vars = "ID")
ind_R <- melt(unique(clean[, .(ID= R, ind_R= median_R)]), id.vars = "ID")
dat <- rbindlist(list(classes, counts, ind_L, ind_R), fill= T)
# marginals
mat <- as.matrix(dcast(clean, L~R, value.var= "diff"), 1)
while(sum(is.na(mat))>0.05*length(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[, -which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}
mar.L <- as.data.table(apply(mat, 1, function(x) mean(x, na.rm= T)), keep.rownames = T)
setnames(mar.L, c("ID", "value"))
mar.L[, variable:= "marginals_L"]
mar.R <- as.data.table(apply(mat, 2, function(x) mean(x, na.rm= T)), keep.rownames = T)
setnames(mar.R, c("ID", "value"))
mar.R[, variable:= "marginals_R"]
dat <- rbind(dat, 
             mar.L,
             mar.R,
             fill= T)
# PCA
imp.L <- apply(mat, 2, function(x) ifelse(is.na(x), median(x, na.rm= T), x))
pca.L <- prcomp(scale(imp.L))
pca.x.L <- as.data.table(pca.L$x, keep.rownames= "ID")
imp.R <- apply(mat, 1, function(x) ifelse(is.na(x), median(x, na.rm= T), x))
pca.R <- prcomp(scale(imp.R))
pca.x.R <- as.data.table(pca.R$x, keep.rownames= "ID")
dat[pca.x.L, c("PC1.L", "PC2.L"):= .(PC1, PC2), on= "ID"]
dat[pca.x.R, c("PC1.R", "PC2.R"):= .(PC1, PC2), on= "ID"]
saveRDS(list(pca.L= pca.L, 
             pca.R= pca.R), 
        "Rdata/vllib002_PCA_LR.rds")

#-----------------------------------------------#
# Plots
#-----------------------------------------------#
Cc_counts <- circlize::colorRamp2(seq(0, 4, length.out= 5),
                                  colorRampPalette(c("cornflowerblue", "gold", "tomato"))(5))
Cc_indAct <- circlize::colorRamp2(seq(0, 8, length.out= 5),
                                  colorRampPalette(c("cornflowerblue", "gold", "tomato"))(5))
marg_palette <- vl_palette_blueWhiteRed(11)
marg_palette[6] <- "grey"
Cc_marg <- circlize::colorRamp2(seq(-1.5, 1.5, length.out= 5),
                                colorRampPalette(marg_palette)(5))

pdf("pdf/draft/PCA_LR_vllib002.pdf",
    width= 9,
    height = 8.5)
par(mfrow= c(3, 3),
    mar= c(3,3,1,1),
    mgp= c(2,0.35,0),
    tcl= -0.2,
    las= 1)
dat[, {
  if(variable=="class")
    .cols <- col else if(grepl("marginals", variable))
    {
      .cols <- Cc_marg(value)
      .breaks <- attr(Cc_marg, "breaks")
      Cv <- Cc_marg(.breaks)
    }else if(grepl("motif$", variable))
    {
      .cols <- Cc_counts(value)
      .breaks <- attr(Cc_counts, "breaks")
      Cv <- Cc_counts(.breaks)
    }else if(grepl("^ind_", variable))
    {
      .cols <- Cc_counts(value)
      .breaks <- attr(Cc_counts, "breaks")
      Cv <- Cc_counts(.breaks)
    }
  
  plot(PC1.L, 
       PC1.R, 
       col= adjustcolor(.cols, 0.8),
       pch= 16,
       main= variable)
  if(variable=="class")
    legend("bottomleft",
           cex= 0.9,
           bty= "n",
           legend= unique(group),
           col= unique(col),
           pch= 16,
           xpd= T) else
           {
             vl_heatkey(.breaks,
                        Cv,
                        left= par("usr")[1]+strwidth("M", cex= 0.5),
                        top= par("usr")[3]+strheight("M")*5,
                        height = strheight("M", cex= 4),
                        main = "Counts",
                        main.cex = 0.9)
           }
  print("")
}, keyby= .(variable)]
dev.off()
