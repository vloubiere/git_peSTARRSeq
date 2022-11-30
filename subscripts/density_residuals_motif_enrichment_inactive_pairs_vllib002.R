setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#---------------------------------#
# Import data
#---------------------------------#
dat <- readRDS("db/FC_tables_DESeq2/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_DESeq2_final_oe.rds")
dat <- dat[actClassL=="inactive" & actClassR=="inactive"]
dat[, meanActL:= mean(log2FoldChange), L]
dat[, meanActR:= mean(log2FoldChange), R]

#---------------------------------#
# Residuals density
#---------------------------------#
pl <- rbind(dat[, .(ID= L, meanAct= meanActL, side= "5' candidates")],
            dat[, .(ID= R, meanAct= meanActR, side= "3' candidates")])
pl <- unique(pl)
pl[, gp:= cut(meanAct,
              c(-Inf, 0, Inf),
              c("Inactive", "Active"),
              include.lowest=T), side]

#---------------------------------#
# Symmetry between left and right
#---------------------------------#
sym <- merge(pl[side== "3' candidates"], 
             pl[side== "5' candidates"], 
             by= "ID", 
             suffixes= c("3'", "5'"))
sym[, gp:= fcase(`gp3'`=="Active" & `gp5'`=="Active", "Active",
                 default = "Inactive")]
sym[, gp:= factor(gp, c("Active", "Inactive"))]
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
# enr <- vl_motif_enrich(counts = counts[gp=="Active", ..sel],
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
mot <- dat[countsL>=3 & countsR>=3]
noMot <- dat[countsL==0 & countsR==0]
noMot <- mot[, {
  idx <- seq(.N)
  res <- .SD[, {
    noMot[, dist:= sqrt((cL-indL)^2+(cR-indR)^2)]
    noMot <- noMot[order(dist)]
    .c <- noMot[1]
    noMot <- noMot[-1]
    .c
  }, .(cL= indL, cR= indR, idx)]
}]
res <- list(noMot$indL,
            mot$indL,
            noMot$indR,
            mot$indR,
            noMot$log2FoldChange,
            mot$log2FoldChange)

#-------------------------------------------------------------#
# PLOT
#-------------------------------------------------------------#
pdf("pdf/draft/density_residuals_motif_enrichment_inactive_pairs_vllib002.pdf", 
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
    .d <- density(meanAct)
    x <- .d$x
    y <- .d$y
    plot(NA, 
         xlim= range(x),
         ylim= range(y),
         ylab= "Density",
         xlab= "mean Obs./Exp. (log2)",
         frame= F,
         main= side)
    polygon(x, y, col= "grey95")
    sup_lim <- min(meanAct[gp=="Active"])
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
  plot(`meanAct3'`,
       `meanAct5'`,
       col= col,
       frame= F,
       pch= 16,
       xlab= "3' candidates mean residuals",
       ylab= "5' candidates mean residuals") 
},]
abline(0, 1, lty= 2)
legend("topleft", 
       legend= paste0("PCC= ", 
                      round(cor.test(sym$`meanAct3'`,
                                     sym$`meanAct5'`)$estimate, 2)),
       bty= "n")

# Twist motif
par(mar= c(6,4,5,6))
xpos <- seq(1, 5, length.out= 3)
vl_boxplot(res,
           at= rep(xpos, each= 2)+c(-0.35, 0.35),
           tilt.names= T, 
           compute_pval= list(c(1,2), c(3,4), c(5,6)),
           main= paste0("Ebox/CATATG/twi", " (", 
                        length(unique(mot$L)), " x 5'; ",
                        length(unique(mot$R)), " x 3'; ",
                        nrow(mot), " pairs)"),
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

# Enrichment
par(mar= c(5,27,2,6),
    mgp= c(3,1,0),
    mfrow= c(1,1),
    las= 2)
padj <- 5e-2
if(nrow(coll[padj<0.05 & log2OR>0]))
{
  epl <- plot(coll, 
              padj_cutoff= padj,
              top_enrich= 17, 
              cex.balloons= 1.5,
              order= "log2OR")
  vl_add_motifs(epl)
}else{
  plot.new()
  text(0.5, 0.5, "No enrichment found (padj<0.05)")
}
dev.off()

saveRDS(unique(epl[, .(name, variable)]),
        "Rdata/top_enrich_motifs_residuals_density_inactive_pairs.rds")