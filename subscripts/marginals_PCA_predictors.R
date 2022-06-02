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
# Sample for CV
set.seed(1)
sel_L <- clean[, 1-.N/clean[,.N], L][, sample(L, round(.N/10), prob = V1)]
set.seed(1)
sel_R <- clean[, 1-.N/clean[,.N], R][, sample(R, round(.N/10), prob = V1)]
clean[, set:= ifelse(L %in% sel_L | R %in% sel_R, "test", "train")]
# Train linear model
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= clean[set=="train"])
# Predict and CV
rsq <- vl_model_eval(observed = clean[set=="test", log2FoldChange], 
                     predicted = predict(model, new= clean[set=="test"]))
clean[, diff:= log2FoldChange-predict(model, new= clean)]

#-----------------------------------------------#
# Format data
#-----------------------------------------------#
# Counts and classes
sub <- feat[ID %in% unique(unlist(clean[, .(L, R)]))]
classes <- sub[, .(ID, group, variable= "class", col)]
counts <- melt(sub[, .(ID, group, Dref_motif, kay_Jra_1_motif, grn_4_motif, grn_24_motif)], id.vars= c("ID", "group"))
counts[, value:= log2(value+1)]
ind_L <- melt(unique(clean[, .(ID= L, ind_L= median_L)]), id.vars = "ID")
ind_R <- melt(unique(clean[, .(ID= R, ind_R= median_R)]), id.vars = "ID")
res <- rbindlist(list(classes, counts, ind_L, ind_R), fill= T)
# Marginals
res[, c("marginals.x", "marginals.y"):= {
  .cL <- clean[.BY, diff, on= "L==ID"]
  .cR <- clean[.BY, diff, on= "R==ID"]
  .(ifelse(length(.cL)>100, mean(.cL), as.numeric(NA)),
    ifelse(length(.cR)>100, mean(.cR), as.numeric(NA)))
}, ID]
# PCA
mat <- as.matrix(dcast(clean, L~R, value.var= "diff"), 1)
while(sum(is.na(mat))>0.05*length(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[, -which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}
imp.L <- apply(mat, 2, function(x) ifelse(is.na(x), median(x, na.rm= T), x))
pca.L <- prcomp(scale(imp.L))
pca.x.L <- as.data.table(pca.L$x, keep.rownames= "ID")[, .(ID, x= PC1, y= PC2)]
imp.R <- apply(mat, 1, function(x) ifelse(is.na(x), median(x, na.rm= T), x))
pca.R <- prcomp(scale(imp.R))
pca.x.R <- as.data.table(pca.R$x, keep.rownames= "ID")[, .(ID, x= PC1, y= PC2)]
res[pca.x.L, c("pca.L.x", "pca.L.y"):= .(i.x, i.y), on= "ID"]
res[pca.x.R, c("pca.R.x", "pca.R.y"):= .(i.x, i.y), on= "ID"]
# melt
dat <- melt(res, 
            measure.vars = patterns("x"= ".x$", "y"= ".y$"),
            variable.name = "coor")
dat[, coor:= c("marginals", "pca.L", "pca.R")[coor]]

#-----------------------------------------------#
# Plots
#-----------------------------------------------#
.palette <- function(x) rev(viridis::viridis(x))

pdf("test/spearman_residuals_marginals.pdf",
    width= 21,
    height = 8.5)
par(mfrow= c(3, 7),
    mar= c(3,3,1,1),
    mgp= c(2,0.35,0),
    tcl= -0.2,
    las= 1)
dat[, {
  if(variable=="class")
    .cols <- col else
    {
      breaks <- seq(0, quantile(value, 0.99), length.out= 10)
      Cv <- .palette(length(breaks))
      Cc <- circlize::colorRamp2(breaks, Cv)
      .cols <- Cc(value)
    }
  plot(x, 
       y, 
       col= adjustcolor(.cols, 0.8),
       pch= 16,
       main= paste(coor, variable))
  if(variable=="class")
    legend("topleft",
           cex= 0.9,
           bty= "n",
           legend= unique(group),
           col= unique(col),
           pch= 16,
           xpd= T) else
           {
             vl_heatkey(breaks,
                        Cv,
                        left= par("usr")[1]+strwidth("M", cex= 0.5),
                        top= par("usr")[4]-strheight("M", cex= 2),
                        height = strheight("M", cex= 4),
                        main = "Counts",
                        main.cex = 0.9)
           }
  print("")
}, keyby= .(coor, variable)]
dev.off()
