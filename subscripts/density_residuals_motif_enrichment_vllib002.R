setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#---------------------------------#
# Import data
#---------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- unique(dat[actClassL!= "inactive", .(L, indL)])[indL<=quantile(indL, 0.9), L]
QR <- unique(dat[actClassR!= "inactive", .(R, indR)])[indR<=quantile(indR, 0.9), R]
dat <- dat[L %in% QL & R %in% QR]

#---------------------------------#
# Residuals density
#---------------------------------#
pl <- rbind(dat[, .(ID= L, mean_residuals= meanResidualsL, side= "5' candidates")],
            dat[, .(ID= R, mean_residuals= meanResidualsR, side= "3' candidates")])
pl <- unique(pl)
pl[, gp:= cut(mean_residuals, 
              quantile(mean_residuals, c(0, 0.1, 0.9, 1)),
              c("Weaker", "Medium", "Stronger"),
              include.lowest=T), side]
pl[, gp:= factor(gp, c("Stronger", "Medium", "Weaker"))]

#---------------------------------#
# Symmetry residuals
#---------------------------------#
sym <- merge(pl[side== "3' candidates"], pl[side== "5' candidates"], by= "ID", suffixes= c("3'", "5'"))

#---------------------------------#
# Motif enrichment
#---------------------------------#
# Define control set
set.seed(1)
ctl <- data.table(side= "None",
                  gp= "ctl",
                  sequence= vl_getSequence(
                    vl_random_regions_BSgenome(genome = "dm3", 
                                               n= 1000,
                                               width = 249),
                    genome= "dm3"))
sel <- vl_Dmel_motifs_DB_full[!is.na(FBgn), motif_ID]
# Counts
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
seq <- pl[, !"mean_residuals"]
seq[lib, sequence:= i.enh_sequence, on= "ID==ID_vl"]
seq$ID <- NULL
# seq <- rbind(seq[gp %in% c("Weaker", "Stronger")], ctl)
seq <- rbind(seq, ctl)
counts <- cbind(seq, vl_motif_counts(seq$sequence, sel= sel))
# Enrichment
enr <- vl_motif_cl_enrich(split(counts[, ..sel], 
                                counts[, .(gp, side)], drop = T), 
                          control_cl = "ctl.None")
setorderv(enr, "padj")
enr <- enr[variable %in% enr[, variable[1], name]$V1] # Select top enrichment

#-------------------------------------------------------------#
# PLOT
#-------------------------------------------------------------#
pdf("pdf/draft/density_residuals_motif_enrichment_vllib002.pdf", 8, height = 8.5)
par(mar= c(10,5,9,1),
    mgp= c(1.75,0.35,0),
    las= 1,
    tcl= -0.2,
    mfrow= c(2,2),
    oma= c(0,0,0,5))
pl[, {
  .c <- unique(.SD)
  .c[, {
    .d <- density(mean_residuals)
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
    sub_lim <- max(mean_residuals[gp=="Weaker"])
    polygon(c(x[x<=sub_lim], sub_lim),
            c(y[x<=sub_lim], 0),
            col= "cornflowerblue")
    sup_lim <- min(mean_residuals[gp=="Stronger"])
    polygon(c(sup_lim, x[x>=sup_lim]),
            c(0, y[x>=sup_lim]),
            col= "tomato")
  }, ]
}, side]
legend(par("usr")[2]-diff(grconvertX(c(0, 1), "line", "user")),
       par("usr")[4],
       legend= rev(levels(pl$gp)),
       fill= c("cornflowerblue", "grey95", "tomato"),
       bty= "n", 
       xpd= NA, 
       x.intersp = 0.2,
       cex= 0.8)

# Symmetry
par(mar= c(6,5,5,1),
    mgp= c(2,0.35,0))
plot(sym$`mean_residuals3'`,
     sym$`mean_residuals5'`,
     col= "lightgrey",
     frame= F,
     pch= 16,
     xlab= "3' candidates mean residuals",
     ylab= "5' candidates mean residuals")
abline(0, 1, lty= 2)
legend("topleft", 
       legend= paste0("PCC= ", 
                      round(cor.test(sym$`mean_residuals3'`,
                                     sym$`mean_residuals5'`)$estimate, 2)),
       bty= "n")
plot.new()

# Enrichment
par(mar= c(14,23,8,6),
    mgp= c(3,1,0),
    mfrow= c(1,1),
    las= 2)
epl <- plot(enr, 
            padj_cutoff= 0.05,
            top_enrich= 7, 
            cex.balloons= 0.9)
vl_add_motifs(epl)
dev.off()

saveRDS(unique(epl[, .(name, variable)]),
        "Rdata/top_enrich_motifs_residuals_density.rds")
