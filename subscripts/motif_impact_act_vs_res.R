setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

dat <- readRDS("db/linear_models/FC_lasso_DSCP_large_WT_residuals_predictions.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
dat$ctlL <- dat$ctlR <- NULL

# Select motifs
ID <- readRDS("db/motif_counts/lib8_motifs_IDs.rds")
mot <- readRDS("db/motif_counts/twist008_motif_counts.rds")
setnames(mot,
         names(mot)[-1],
         as.character(ID$cluster))
dat <- merge(dat, mot, by.x= "L", by.y= "ID", all.x= T, allow.cartesian = T)
dat <- merge(dat, mot, by.x= "R", by.y= "ID", all.x= T, allow.cartesian = T, suffixes= c("__L", "__R"))

# Compute residuals/activity for each motif combination and t.test pval ----
cmb <- CJ(grep("__L$", names(dat), value = T),
          grep("__R$", names(dat), value = T))
res <- cmb[, {
  if(.GRP %% 100==0)
    print(.GRP)
  # Check motifs
  mot <- dat[[V1]]>0 & dat[[V2]]>0 # Expected motifs on both sides
  act_mot <- dat[(mot), log2FoldChange]
  res_mot <- dat[(mot), residuals]
  noMot <- dat[[V1]]==0 & dat[[V2]]==0 # No motifs
  act_nomot <- dat[(noMot), log2FoldChange]
  res_nomot <- dat[(noMot), residuals]
  # Return values
  .(act_mot= mean(act_mot),
    act_nomot= mean(act_nomot),
    act_pval= t.test(act_mot, act_nomot)$p.value,
    res_mot= mean(res_mot),
    res_nomot= mean(res_nomot),
    res_pval= t.test(res_mot, res_nomot)$p.value)
}, .(V1, V2)]
# FDR ----
res[, act_padj:= p.adjust(act_pval, "fdr")]
res[, res_padj:= p.adjust(res_pval, "fdr")]
plot(density(-log10(res$res_padj)),
     xlab= "-log10(padj)",
     col= "red",
     frame= F)
lines(density(-log10(res$act_padj)))
legend("topright",
       col= c("red", "black"),
       legend= c("Residuals", "Activity"),
       lwd= 1,
       bty= "n")
# Simplify names ----
res[, V1:= gsub("/", ".", V1)]
res[, V2:= gsub("/", ".", V2)]
# Compute difference between pairs with motifs and without motifs ----
res[, `Mean activity [motif / no motif] (log2)`:= act_mot-act_nomot]
res[, `Mean residuals [motif / no motif] (log2)`:= res_mot-res_nomot]
# Add colors for pairs of interest ----
res[, sel:= V1 %in% c("Ebox.CATATG.twi__L", "Trl.1__L") & V2 %in% c("Ebox.CATATG.twi__R", "Trl.1__R")
    | (V1 %in% c("DRE.1__L", "DRE.2__L") & V2 %in% c("DRE.1__R", "DRE.2__R"))
    | (V1 %in% c("AP1.1__L", "GATA.1__L") & V2 %in% c("AP1.1__R", "GATA.1__R"))]
res[, col:= fcase(V1 %in% c("Ebox.CATATG.twi__L", "Trl.1__L") & V2 %in% c("Ebox.CATATG.twi__R", "Trl.1__R"), "tomato",
                  V1 %in% c("DRE.1__L", "DRE.2__L") & V2 %in% c("DRE.1__R", "DRE.2__R"), "cornflowerblue",
                  V1 %in% c("AP1.1__L", "GATA.1__L") & V2 %in% c("AP1.1__R", "GATA.1__R"), "limegreen",
                  default = "lightgrey")]
res[, label:= paste0(gsub("__L$", "", V1), " / ", gsub("__R$", "", V2))]
res[, label:= gsub("AP1.1", "AP-1", label)]
res[, label:= gsub("GATA.1", "GATA", label)]
res[, label:= gsub("Ebox.CATATG.twi", "Twist", label)]
res[, label:= gsub("DRE.1|DRE.2", "Dref", label)]
res[, label:= gsub("Trl.1", "Trl", label)]
res[, label:= gsub("SREBP.1", "SREBP", label)]
res[, label:= gsub("CREB.ATF.2", "CREB", label)]
setorderv(res, "sel")

# Heatmap to compare symmetry ----
res[, V1:= gsub("__L$", "", V1)]
res[, V2:= gsub("__R$", "", V2)]
mat <- dcast(res,
             V1~V2,
             value.var = "Mean residuals [motif / no motif] (log2)")
mat <- as.matrix(mat, 1)

# Add examples to boxplot homotypic motif pairs ----
pl <- function(left, right, x, col, xadj= .1, yadj= 0.05)
{
  .c <- res[V1 %in% left & V2 %in% right]
  y <- .c[, `Mean residuals [motif / no motif] (log2)`]
  x <- rep(x, length(y))
  segments(x, y, x+xadj, y+yadj)
  text(x+xadj,
       y+yadj,
       .c$label,
       pos= 4,
       xpd= T,
       cex= 5/12,
       offset= .1)
  points(x,
         y,
         bg= col,
         pch= 21,
         cex= .5)
}

# Plot ----
pdf("pdf/draft/motif_impact_act_vs_res.pdf", 
    width = 3,
    height = 3)
par(mai= c(0.75,0.75,0.75,0.75), 
    mgp= c(1, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
plot.new()
plot.window(xlim= c(-4.5,4.5),
            ylim= c(0,.5))
rect(par("usr")[1],
     0,
     -2,
     .5,
     col= adjustcolor("cornflowerblue", .3),
     border= NA)
text(mean(c(par("usr")[1], -2)),
     .5,
     pos= 1,
     labels = "Weaker",
     cex= 7/12)
rect(2,
     0,
     par("usr")[2],
     .5,
     col= adjustcolor("tomato", .3),
     border= NA)
text(mean(c(2, par("usr")[2])),
     .5,
     pos= 1,
     labels = "Stronger",
     cex= 7/12)
par(lwd= 0.5)
hist(dat$residuals,
     xlim= c(-4.5,4.5),
     freq = F,
     add= T)
axis(1, padj = -1.25)
axis(2)
title(xlab= "Residuals (log2)")
title(ylab= "Density")
res[, {
  vl_repelScatterplot(`Mean activity [motif / no motif] (log2)`,
                      `Mean residuals [motif / no motif] (log2)`,
                      labels = label,
                      label.cex = 5/12,
                      label.sel = sel,
                      col= adjustcolor(col, .4),
                      pch= 16,
                      ylim= c(-1, 1),
                      cex= .5)
}]
# Boxplot homotypic vs heterotypic
vl_par(mai = c(.7, 1, .7, 1))
vl_boxplot(res[V1==V2, `Mean residuals [motif / no motif] (log2)`],
           res[V1!=V2, `Mean residuals [motif / no motif] (log2)`],
           compute.pval = list(c(1,2)),
           names= c("Homotypic pairs",
                    "Heterotypic pairs"),
           tilt.names = T,
           ylab= "Mean residuals\n[motif / no motif] (log2)")
pl("AP1.1", "AP1.1", 1, "tomato")
pl("AP1.1",
   c("CREB.ATF.2", "GATA.1", "SREBP.1", "Ebox.CATATG.twi"),
   2,
   "limegreen",
   yadj= c(.05, .0, .05, .05))
dev.off()