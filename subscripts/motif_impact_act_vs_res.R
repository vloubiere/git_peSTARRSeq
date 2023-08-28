setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
dat$ctlL <- dat$ctlR <- NULL

# Select motifs
ID <- readRDS("db/motif_counts/motifs_IDs.rds")
mot <- readRDS("db/motif_counts/twist008_motif_counts.rds")
setnames(mot,
         names(mot)[-1],
         as.character(ID$cluster))
dat <- merge(dat, mot, by.x= "L", by.y= "ID", all.x= T, allow.cartesian = T)
dat <- merge(dat, mot, by.x= "R", by.y= "ID", all.x= T, allow.cartesian = T, suffixes= c("__L", "__R"))

# Compute residuals for each motif combination
cmb <- CJ(grep("__L$", names(dat), value = T),
          grep("__R$", names(dat), value = T))
res <- cmb[, {
  if(.GRP %% 100==0)
    print(.GRP)
  .(act_mot= mean(dat$log2FoldChange[dat[[V1]]>0 & dat[[V2]]>0]),
    act_nomot= mean(dat$log2FoldChange[dat[[V1]]==0 & dat[[V2]]==0]),
    res_mot= mean(dat$residuals[dat[[V1]]>0 & dat[[V2]]>0]),
    res_nomot= mean(dat$residuals[dat[[V1]]==0 & dat[[V2]]==0]))
}, .(V1, V2)]
res[, `Mean activity [motif>0/motif=0] (log2)`:= act_mot-act_nomot]
res[, `Mean residuals [motif>0/motif=0] (log2)`:= res_mot-res_nomot]
res[, V1:= gsub("__L$", "", V1)]
res[, V2:= gsub("__R$", "", V2)]
res[, col:= fcase(V1 %in% c("Trl/1", "Ebox/CATATG/twi"), "tomato",
                  V1 %in% c("DRE/1"), "cornflowerblue",
                  default = "lightgrey")]
setorderv(res, "res_mot", -1)
res[, col1:= fcase(.I %in% c(1,3,6,19), "tomato",
                   .I==.N, "cornflowerblue",
                   default = "lightgrey")]
mat <- dcast(res,
             V1~V2,
             value.var = "Mean residuals [motif>0/motif=0] (log2)")
mat <- as.matrix(mat, 1)

# Plots
pdf1 <- tempfile(fileext = ".pdf")
pdf2 <- tempfile(fileext = ".pdf")
pdf3 <- tempfile(fileext = ".pdf")

## Heatmap 
pdf(pdf1, 4, 4)
par(mar= c(15,11.5,2,7.5),
    las= 2,
    cex= 0.3,
    lwd= 0.5)
vl_heatmap(mat,
           breaks = c(-0.25, 0, 0.25),
           legend_title = "Mean\nresiduals\n[motif>0/motif=0]\n(log2)\n")
dev.off()

## Act vs res
pdf(pdf2, 3.5, 3.5)
par(mar= c(4,4,3,3),
    mgp= c(2.5,0.5,0),
    las= 1,
    lwd= 0.5,
    cex= 0.8)
res[V1==V2, {
  vl_repelScatterplot(`Mean activity [motif>0/motif=0] (log2)`,
                      `Mean residuals [motif>0/motif=0] (log2)`,
                      pch= 16,
                      col= col,
                      labels = V1,
                      rect.col = adjustcolor(col, 0.5),
                      label.cex = ifelse(col=="lightgrey", .3, .7),
                      point_padding_x = 0.01, 
                      point_padding_y = 0.01)
}]
dev.off()

# Sorted scatterplot
pdf(pdf3, 7, 3.5)
par(las= 1,
    tcl= -0.2,
    mgp= c(2, 0.5, 0),
    bty= "n")
layout(matrix(1:3, nrow= 1), widths = c(2/5, 1/5, 2/5))
dat[, {
  smoothScatter(predicted,
                residuals,
                xlab= "Fitted values (log2)",
                ylab= "Residuals (log2)",
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                nrpoints = 0)
  abline(h= 0, 
         lwd= 0.5)
}]
dat[, {
  vl_boxplot(residuals,
             ylab= "Residuals (log2)",
             violin= T,
             col= "white",
             viocol= "lightgrey")
  abline(h= 0,
         lwd= 0.5)
}]
res[, {
  vl_repelScatterplot(res_mot,
                      pch= 16,
                      col= adjustcolor(col1, .5),
                      rect.col= adjustcolor(col1, .5),
                      labels = paste0(V1, " - ", V2),
                      label.sel = col1!="lightgrey",
                      xaxt= "n",
                      ylab= "Mean residuals (log2)",
                      xlab= NA)
  abline(h= 0, lwd= 0.5)
}]
dev.off()

pdftools::pdf_combine(c(pdf1, pdf2, pdf3),
                      "pdf/draft/motif_impact_act_vs_res.pdf")