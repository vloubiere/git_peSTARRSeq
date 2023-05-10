setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#---------------------------------#
# Import data
#---------------------------------#
dat <- readRDS("db/linear_models/FC_vllib030_with_predictions.rds")
dat <- dat[!grepl("control", L) & !grepl("control", R)]

#---------------------------------#
# Residuals density
#---------------------------------#
pl <- rbind(dat[, .(ID= L, mean_residuals= meanResidualsL, side= "5' candidates")],
            dat[, .(ID= R, mean_residuals= meanResidualsR, side= "3' candidates")])
pl <- unique(pl)
pl[, gp:= cut(mean_residuals,
              quantile(mean_residuals, c(0, 0.25, 0.75, 1)),
              c("Sub-efficient", "Synergistic", "Super-efficient"),
              include.lowest=T), side]
pl[, gp:= factor(gp, c("Super-efficient", "Synergistic", "Sub-efficient"))]

#---------------------------------#
# Symmetry residuals
#---------------------------------#
sym <- merge(pl[side=="3' candidates"],
             pl[side=="5' candidates"],
             by= "ID",
             suffixes= c("5'", "3'"))
sym[, gp:= `gp5'`]
sym[`gp5'`!=`gp3'`, gp:= "Synergistic"]
sym[, col:= c("tomato", "lightgrey", "cornflowerblue")[gp]]

#---------------------------------#
# Motif enrichment
#---------------------------------#
# Counts
counts <- readRDS("db/motif_counts/twist015_motif_counts_low_stringency_no_collapsing.rds")
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
pdf("pdf/draft/density_residuals_motif_enrichment_vllib030.pdf",
    width= 8,
    height = 9)
par(mar= c(6,22,8,22),
    mgp= c(1.75,0.35,0),
    las= 1,
    tcl= -0.2,
    mfrow= c(3,1))
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
    legend(par("usr")[2]-diff(grconvertX(c(0, 1), "line", "user")),
           par("usr")[4],
           legend= rev(levels(gp)),
           fill= c("cornflowerblue", "grey95", "tomato"),
           bty= "n",
           xpd= NA)
  }, ]
}, side]

# Symmetry
par(mar= c(4,22,3,22),
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
par(mar= c(11,27,11,10),
    mgp= c(3,1,0),
    mfrow= c(1,1),
    las= 2)
epl <- plot(coll, 
            padj_cutoff= 0.05,
            cex.balloons= 0.9, 
            order= "log2OR", 
            top_enrich= 8)
vl_add_motifs(epl)
dev.off()

# saveRDS(unique(epl[, .(name, variable)]),
#         "Rdata/top_enrich_motifs_residuals_density.rds")

file.show("pdf/draft/density_residuals_motif_enrichment_vllib030.pdf")
