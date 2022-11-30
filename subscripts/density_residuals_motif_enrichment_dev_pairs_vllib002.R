setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#---------------------------------#
# Import data
#---------------------------------#
dat <- readRDS("db/linear_models/FC_dev_pairs_vllib002_with_predictions.rds")

#---------------------------------#
# Residuals density
#---------------------------------#
pl <- rbind(dat[, .(ID= L, mean_residuals= meanResidualsL, side= "5' candidates")],
            dat[, .(ID= R, mean_residuals= meanResidualsR, side= "3' candidates")])
pl <- unique(pl)
pl[, gp:= cut(mean_residuals,
              quantile(mean_residuals, c(0, 0.33, 0.66, 1)),
              c("Sub-efficient", "Synergistic", "Super-efficient"),
              include.lowest=T), side]
pl[, gp:= factor(gp, c("Super-efficient", "Synergistic", "Sub-efficient"))]

#---------------------------------#
# Symmetry residuals
#---------------------------------#
sym <- merge(pl[side== "3' candidates"], 
             pl[side== "5' candidates"],
             by= "ID", 
             suffixes= c("3'", "5'"))
sym[, gp:= fcase(`gp3'`=="Sub-efficient" & `gp5'`=="Sub-efficient", "Sub-efficient",
                 `gp3'`=="Super-efficient" & `gp5'`=="Super-efficient", "Super-efficient",
                 default = "Synergistic")]
sym[, gp:= factor(gp, c("Super-efficient", "Synergistic", "Sub-efficient"))]
sym[, col:= c("tomato", "lightgrey", "cornflowerblue")[gp]]

#---------------------------------#
# Motif enrichment
#---------------------------------#
# Counts
counts <- readRDS("db/motif_counts/twist008_motif_counts_low_stringency_no_collapsing.rds")
sel <- names(counts)[-1]
counts <- counts[sym$ID, on= "ID"]
ctls <- readRDS("db/motif_counts/random_controls_1000_low_stringency_no_collapsing.rds")
# Enrichment
enr <- vl_motif_cl_enrich(c(split(counts[, ..sel], sym$gp),
                            list(ctl= ctls[, ..sel])), 
                          control_cl = "ctl")
setorderv(enr, "padj")
coll <- enr[variable %in% enr[, variable[1], name]$V1] # Select top enrichment

#-------------------------------------------------------------#
# PLOT
#-------------------------------------------------------------#
pdf("pdf/draft/density_residuals_motif_enrichment_vllib002.pdf",
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
    sub_lim <- max(mean_residuals[gp=="Sub-efficient"])
    polygon(c(x[x<=sub_lim], sub_lim),
            c(y[x<=sub_lim], 0),
            col= "cornflowerblue")
    sup_lim <- min(mean_residuals[gp=="Super-efficient"])
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
sym[, {
  plot(`mean_residuals3'`,
       `mean_residuals5'`,
       col= col,
       frame= F,
       pch= 16,
       xlab= "3' candidates mean residuals",
       ylab= "5' candidates mean residuals") 
},]
abline(0, 1, lty= 2)
legend("topleft", 
       legend= paste0("PCC= ", 
                      round(cor.test(sym$`mean_residuals3'`,
                                     sym$`mean_residuals5'`)$estimate, 2)),
       bty= "n")

# Enrichment
plot.new()
par(mar= c(8,25,2,6),
    mgp= c(3,1,0),
    mfrow= c(1,1),
    las= 2)
epl <- plot(coll, 
            padj_cutoff= 0.001,
            top_enrich= 9,
            cex.balloons= 0.9, 
            order= "log2OR")
vl_add_motifs(epl)
dev.off()

saveRDS(unique(epl[, .(name, variable)]),
        "Rdata/top_enrich_motifs_residuals_density.rds")

