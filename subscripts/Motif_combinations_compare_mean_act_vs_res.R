setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

# Import mean act/residuals per motif combination ----
dat <- readRDS("db/motif_effect/motif_impact_act_res.rds")

# Add names ----
names <- readRDS("db/motif_counts/lib8_motifs_enrichments.rds")
cols <- c("V1", "V2")
dat[, (cols):= lapply(.SD, function(x) gsub("__L$|__R$", "", x)), .SDcols= cols]

# Add colors for pairs of interest ----
sel <- data.table(name= c("GATA", "AP-1", "Twist", "Trl",  "Dref", "Dref"),
                  motif_ID= c("cisbp__M4320", "cisbp__M6317",
                              "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715",
                              "homer__CTCTCTCTCY_GAGA-repeat",
                              "flyfactorsurvey__Dref_FlyReg_FBgn0015664",
                              "homer__AVYTATCGATAD_DREF"),
                  col= rep(c("limegreen", "tomato", "cornflowerblue"), each= 2))
sel <- merge(sel, sel, by= "col")
sel[, label:= paste0(name.x, "/", name.y)]
dat <- merge(dat, sel, by.x= c("V1", "V2"), by.y= c("motif_ID.x", "motif_ID.y"), all.x= T)
dat[is.na(col), col:= "lightgrey"]
dat <- dat[order(!is.na(label))]

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
dat[, {
  vl_repelScatterplot(`Mean activity [motif / no motif] (log2)`,
                      `Mean residuals [motif / no motif] (log2)`,
                      labels = label,
                      label.cex = 5/12,
                      label.sel = !is.na(label),
                      col= adjustcolor(col, .4),
                      point_size = .1,
                      pch= 16,
                      ylim= c(-1, 1),
                      cex= .5)
}]
dev.off()