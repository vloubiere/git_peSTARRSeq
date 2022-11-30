setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#---------------------------------#
# Import data
#---------------------------------#
dat <- readRDS("db/FC_tables_DESeq2/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_DESeq2_final_oe.rds")
dat <- dat[actClassL=="inactive" & actClassR=="active" | actClassL=="active" & actClassR=="inactive"]
dat[, meanDiffL:= mean(log2FoldChange-indR), L]
dat[, meanDiffR:= mean(log2FoldChange-indL), R]

#---------------------------------#
# Residuals density
#---------------------------------#
pl <- rbind(dat[actClassL=="inactive", .(ID= L, meanDiff= meanDiffL, side= "5' candidates")],
            dat[actClassR=="inactive", .(ID= R, meanDiff= meanDiffR, side= "3' candidates")])
pl <- unique(pl)
pl[, gp:= cut(meanDiff,
              c(-Inf, 0, Inf),
              c("Inactive", "Boosting"),
              include.lowest=T), side]

#---------------------------------#
# Symmetry between left and right
#---------------------------------#
sym <- merge(pl[side== "3' candidates"],
             pl[side== "5' candidates"],
             by= "ID",
             suffixes= c("3'", "5'"))
sym[, gp:= fcase(`gp3'`=="Boosting" & `gp5'`=="Boosting", "Boosting",
                 default = "Inactive")]
sym[, gp:= factor(gp, c("Boosting", "Inactive"))]
sym[, col:= c("tomato", "lightgrey")[gp]]

#---------------------------------#
# Motif enrichment
#---------------------------------#
# Define control set
set.seed(1)
ctl <- data.table(gp= "ctl",
                  sequence= vl_getSequence(
                    vl_random_regions_BSgenome(genome = "dm3", 
                                               n= 1000,
                                               width = 249),
                    genome= "dm3"))
# Retrieve enhancer sequences
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
sym[lib, sequence:= i.enh_sequence, on= "ID==ID_vl"]
# Counts
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
sel <- as.data.table(TF_clusters_PWMs[["metadata"]])[S2_exp>0, motif_name]
seq <- rbind(sym[, .(gp, sequence)], ctl)
counts <- cbind(seq, vl_motif_counts(seq$sequence, sel= sel, genome= "dm3"))
# Enrichment
enr <- vl_motif_cl_enrich(split(counts[, ..sel],
                                counts$gp,
                                drop = T),
                          control_cl = "ctl")
# enr <- vl_motif_enrich(counts = counts[gp=="Boosting", ..sel],
#                        control_counts = counts[gp=="Inactive", ..sel])
setorderv(enr, "padj")
coll <- enr[variable %in% enr[, variable[1], name]$V1] # Select top enrichment

#---------------------------------#
# Twist motif activity (activity matched controls)
#---------------------------------#
twist <- data.table(ID= sym$ID,
                    vl_motif_counts(sym$sequence, 
                                    sel= "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715", 
                                    genome= "dm3", 
                                    collapse_overlapping = T))
dat[twist, countsL:= i.flyfactorsurvey__CG16778_SANGER_5_FBgn0003715, on= "L==ID"]
dat[twist, countsR:= i.flyfactorsurvey__CG16778_SANGER_5_FBgn0003715, on= "R==ID"]
# Left
motL <- dat[countsL>=3 & actClassL=="inactive"]
noMotL <- dat[countsL==0 & actClassL=="inactive"]
noMotL <- motL[, {
  idx <- seq(.N)
  res <- .SD[, {
    noMotL[, dist:= sqrt((cL-indL)^2+(cR-indR)^2)]
    noMotL <- noMotL[order(dist)]
    .c <- noMotL[1]
    noMotL <- noMotL[-1]
    .c
  }, .(cL= indL, cR= indR, idx)]
}]
resL <- list(noMotL$indL,
             motL$indL,
             noMotL$indR,
             motL$indR,
             noMotL$log2FoldChange,
             motL$log2FoldChange)
# Right
motR <- dat[countsR>=3 & actClassR=="inactive"]
noMotR <- dat[countsR==0 & actClassR=="inactive"]
noMotR <- motR[, {
  idx <- seq(.N)
  res <- .SD[, {
    noMotR[, dist:= sqrt((cL-indL)^2+(cR-indR)^2)]
    noMotR <- noMotR[order(dist)]
    .c <- noMotR[1]
    noMotR <- noMotR[-1]
    .c
  }, .(cL= indL, cR= indR, idx)]
}]
resR <- list(noMotR$indL,
             motR$indL,
             noMotR$indR,
             motR$indR,
             noMotR$log2FoldChange,
             motR$log2FoldChange)

#-------------------------------------------------------------#
# PLOT
#-------------------------------------------------------------#
pdf("pdf/draft/density_residuals_motif_enrichment_inactive_active_pairs_vllib002.pdf", 
    width= 8,
    height = 9)
par(mar= c(10,5,9,1),
    mgp= c(1.75,0.35,0),
    las= 1,
    tcl= -0.2,
    mfrow= c(2,2),
    oma= c(0,0,0,5))
pl[, {
  .c <- unique(.SD)
  .c[, {
    .d <- density(meanDiff)
    x <- .d$x
    y <- .d$y
    plot(NA, 
         xlim= range(x),
         ylim= range(y),
         ylab= "Density",
         xlab= "mean Obs./Exp. (log2)",
         frame= F,
         main= side)
    polygon(x, y, col= "lightgrey")
    sup_lim <- min(meanDiff[gp=="Boosting"])
    polygon(c(sup_lim, x[x>=sup_lim]),
            c(0, y[x>=sup_lim]),
            col= "tomato")
  }, ]
}, side]
legend(par("usr")[2]-diff(grconvertX(c(0, 1), "line", "user")),
       par("usr")[4],
       legend= rev(levels(pl$gp)),
       fill= c("tomato", "lightgrey"),
       bty= "n", 
       xpd= NA, 
       x.intersp = 0.2,
       cex= 0.8)

# Symmetry
par(mar= c(6,4,5,1),
    mgp= c(2,0.35,0))
sym[, {
  plot(`meanDiff3'`,
       `meanDiff5'`,
       col= col,
       frame= F,
       pch= 16,
       xlab= "3' candidates mean residuals",
       ylab= "5' candidates mean residuals") 
},]
abline(0, 1, lty= 2)
legend("topleft", 
       legend= paste0("PCC= ", 
                      round(cor.test(sym$`meanDiff3'`,
                                     sym$`meanDiff5'`)$estimate, 2)),
       bty= "n")

# Enrichment
plot.new()
par(mar= c(12,25,12,8),
    mgp= c(3,1,0),
    mfrow= c(1,1),
    las= 2)
padj <- 5e-2
if(nrow(coll[padj<0.05 & log2OR>0]))
{
  epl <- plot(coll, 
              padj_cutoff= padj,
              top_enrich= 20, 
              cex.balloons= 1.5)
  vl_add_motifs(epl)
}else{
  plot.new()
  text(0.5, 0.5, "No enrichment found (padj<0.05)")
}

# Twist motif
par(mar= c(6,7,5,5),
    mfrow= c(2,2),
    las= 1)
for(pos in list(list(mot= motL, res= resL), 
                list(mot= motR, res= resR)))
{
  .c <- pos[["mot"]]
  res <- pos[["res"]]
  xpos <- seq(1, 5, length.out= 3)
  vl_boxplot(res,
             at= rep(xpos, each= 2)+c(-0.35, 0.35),
             tilt.names= T,
             compute_pval= list(c(1,2), c(3,4), c(5,6)),
             main= paste0("Ebox/CATATG/twi", " (", 
                          length(unique(.c$L)), " x 5'; ",
                          length(unique(.c$R)), " x 3'; ",
                          length(res[[4]]), " pairs)"),
             xaxt= "n",
             col= c("lightgrey", "rosybrown1"),
             ylab= "Activity (log2)")
  axis(1, at= xpos, labels = c("5'", "3'", "Pair"))
  legend(par("usr")[2],
         par("usr")[4],
         fill= c("lightgrey", "rosybrown1"),
         legend= c("Controls (no motifs)",
                   paste0("Motif pairs (>= 3)")),
         bty= "n",
         xpd= NA)
}

# Combined boxplot
par(mar= c(6,8,5,6))
xpos <- c(1, 4)
vl_boxplot(noMotL$meanDiffL, 
           motL$meanDiffL,
           noMotR$meanDiffR, 
           motR$meanDiffR, 
           at= rep(xpos, each= 2)+c(-0.35, 0.35),
           compute_pval= list(c(1,2), c(3,4)),
           col= c("lightgrey", "rosybrown1"),
           ylab= c("Pairs' activity - active element's activity (log2)"),
           xaxt= "n",
           main= "Ebox/CATATG/twi")
axis(1, 
     at= xpos, 
     labels = c("5' twist", "3' twist"))
legend(par("usr")[2],
       par("usr")[4],
       fill= c("lightgrey", "rosybrown1"),
       legend= c("Controls (no motifs)",
                 paste0("Motif pairs (>= 3)")),
       bty= "n",
       xpd= NA)
abline(h= 0, lty= 2)
dev.off()

saveRDS(unique(epl[, .(name, variable)]),
        "Rdata/top_enrich_motifs_residuals_density_inactive_active_pairs.rds")