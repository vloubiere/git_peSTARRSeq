setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Train activity based models and compute residuals
#-----------------------------------------------#
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
screen <- readRDS("Rdata/final_results_table.rds")
screen <- screen[vllib=="vllib002" 
                 & class== "enh./enh."]
screen <- feat$add_feature(screen, feat$lib)
screen <- screen[group_L %in% c("hk", "dev", "shared")
                 & group_R %in% c("hk", "dev", "shared")]
# Sample for CV
set.seed(1)
sel_L <- screen[, 1-.N/screen[,.N], L][, sample(L, round(.N/10), prob = V1)]
set.seed(1)
sel_R <- screen[, 1-.N/screen[,.N], R][, sample(R, round(.N/10), prob = V1)]
screen[, set:= ifelse(L %in% sel_L | R %in% sel_R, "test", "train")]
# Train linear model
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= screen[set=="train"])
# Predict and CV
rsq <- vl_model_eval(observed = screen[set=="test", log2FoldChange], 
                     predicted = predict(model, new= screen[set=="test"]))
screen[, diff:= log2FoldChange-predict(model, new= screen)]

#-----------------------------------------------#
# Format data
#-----------------------------------------------#
# Import
dat <- data.table(ID= unique(unlist(screen[, .(L, R)])))
dat[, class:= ]
cl <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")
cl$rows
cl$cols
counts <- readRDS("Rdata/uniq_enh_feat/motifCounts.rds")
counts <- as.data.table(counts)[, .(Dref, kay_Jra_1, grn_4, grn_24)]
setnames(counts, paste0(names(counts), "_motif"))
lib <- readRDS(feat$lib)
lib <- cbind(lib[, .(ID, group, col)], counts)
lib <- lib[group %in% c("hk", "shared", "dev")]



#-----------------------------------------------#
# Format data
#-----------------------------------------------#
dat <- lib[ID %in% unique(unlist(screen[, .(L, R)]))]
# Marginals
dat[, marginal.x:= mean(screen[.BY, diff, on= "L==ID"]), ID]
dat[, marginal.x:= ifelse(screen[.BY, .N, on= "L==ID"]>100, marginal.x, as.numeric(NA)), ID]
dat[, marginal.y:= mean(screen[.BY, diff, on= "R==ID"]), ID]
dat[, marginal.y:= ifelse(screen[.BY, .N, on= "R==ID"]>100, marginal.y, as.numeric(NA)), ID]
# PCA
mat <- as.matrix(dcast(screen, L~R, value.var= "diff"), 1)
while(sum(is.na(mat))>0.05*length(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[, -which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}
imp.L <- apply(mat, 2, function(x) ifelse(is.na(x), median(x, na.rm= T), x))
pca.L <- prcomp(scale(imp.L))
imp.R <- apply(mat, 1, function(x) ifelse(is.na(x), median(x, na.rm= T), x))
pca.R <- prcomp(scale(imp.R))
dat[as.data.table(pca.L$x[, c("PC1", "PC2")], keep.rownames = T), c("PCA_L.x", "PCA_L.y"):= .(PC1, PC2), on= "ID==rn"]
dat[as.data.table(pca.R$x[, c("PC1", "PC2")], keep.rownames = T), c("PCA_R.x", "PCA_R.y"):= .(PC1, PC2), on= "ID==rn"]

#-----------------------------------------------#
# Plotting object
#-----------------------------------------------#
pl <- melt(dat, measure.vars = patterns("x"= ".x$", "y"= ".y$"))
pl[, variable:= c("marginal", "PCA_L", "PCA_R")[variable]]
pl[, col:= factor(col)]
leg <- unique(pl[, .(group, col)])
setkeyv(leg, "group")
pl[, class:= as.numeric(col)]
pl$col <- NULL
pl <- melt(pl, id.vars = c("ID", "group", "x", "y", "variable"), variable.name = "var")
pl[, col:= {
  if(var=="class")
    .c <- cols[value]
  if(grepl("_motif$", var))
  {
    breaks <- seq(0, quantile(value, 0.99))
    .cols <- rev(viridis::viridis(length(breaks)))
    Cc <- circlize::colorRamp2(breaks, .cols)
    .c <- Cc(value)
  }
  .c
}, .(var, variable)]

#-----------------------------------------------#
# Plots
#-----------------------------------------------#
pdf("test/spearman_residuals_marginals.pdf",
    width= 15,
    height = 8.5)
par(mfrow= c(3, 5),
    mar= c(3,3,1,1),
    mgp= c(2,0.35,0),
    tcl= -0.2,
    las= 1)
pl[, {
  plot(x, 
       y, 
       col= adjustcolor(col, 0.8),
       pch= 16,
       main= paste(variable, var))
  if(var=="class")
    legend("topleft",
           cex= 0.9,
           bty= "n",
           legend= leg$group,
           col= adjustcolor(leg$col, 0.8),
           pch= 16,
           xpd= T)
  if(grepl("_motif$", var))
  {
    breaks <- seq(0, quantile(value, 0.99))
    .cols <- rev(viridis::viridis(length(breaks)))
    vl_heatkey(breaks,
               .cols,
               left= par("usr")[1]+strwidth("M", cex= 0.5),
               top= par("usr")[4]-strheight("M", cex= 2),
               height = strheight("M", cex= 4),
               main = "Counts",
               main.cex = 0.9)
  }
  print("")
}, keyby= .(variable, var)]
dev.off()
 
# plot(marg[, .(mar_L, mar_R)],
#      col= adjustcolor(lib$col[match(marg$enh, lib$ID)], 0.7),
#      pch= 16,
#      main= "Marginals") 
# leg()
# abline(0,1)
# abline(1, 1, lty= "11")
# abline(-1, 1, lty= "11")
# for(var in vars)
# {
#   breaks <- seq(0, quantile(counts[ID %in% marg$enh][[var]], 0.99))
#   cols <- rev(viridis::viridis(length(breaks)))
#   Cc <- circlize::colorRamp2(breaks, cols)
#   plot(marg[, .(mar_L, mar_R)],
#        col= Cc(counts[match(marg$enh, counts$ID), var, with= F]),
#        pch= 16,
#        main= "Marginals") 
#   vl_heatkey(breaks, 
#              cols, 
#              left= par("usr")[1]+strwidth("M", cex= 0.5),
#              top= par("usr")[4]-strheight("M", cex= 2), 
#              height = strheight("M", cex= 4), 
#              main = "Counts", 
#              main.cex = 0.9)
# }
# for(side in c("pca.L", "pca.R"))
# {
#   .c <- get(side)$x 
#   plot(.c[, c("PC1", "PC2")],
#        col= adjustcolor(lib$col[match(rownames(.c), lib$ID)], 0.7),
#        pch= 16,
#        main= fcase(side== "pca.L", "5'PCA",
#                    side== "pca.R", "3'PCA"))
#   leg()
#   for(var in vars)
#   {
#     breaks <- seq(0, quantile(counts[ID %in% rownames(.c)][[var]], 0.99))
#     cols <- rev(viridis::viridis(length(breaks)))
#     Cc <- circlize::colorRamp2(breaks, cols)
#     plot(.c[, c("PC1", "PC2")],
#          col= Cc(counts[match(rownames(.c), counts$ID), var, with= F]),
#          pch= 16,
#          main= var)
#     vl_heatkey(breaks, 
#                cols, 
#                left= par("usr")[1]+strwidth("M", cex= 0.5),
#                top= par("usr")[4]-strheight("M", cex= 2), 
#                height = strheight("M", cex= 4), 
#                main = "Counts", 
#                main.cex = 0.9)
#   }
# }
# dev.off()
# 
