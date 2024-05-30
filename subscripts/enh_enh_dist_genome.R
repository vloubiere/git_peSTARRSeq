setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- fread("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/peaks/DSCP_200bp_gw.UMI_cut_merged.peaks.txt",
             sel= c(1,2,5,6,7),
             col.names = c("seqnames", "start", "Enrch.", "Corr_enrch", "p_value"))
dat[, end:= start]
dist <- vl_closestBed(dat, min.dist = 1)

# Plot ----
pdf("pdf/draft/enhancer_enhancer_distance.pdf",
    height = 3,
    width = 2.25)
vl_par(cex.lab= 8/12,
       cex.axis= 7/12)
dist[, {
  vl_boxplot(abs(dist)/1000,
             notch= T,
             ylab = "Enh.-enh. distance (kb)",
             col= "lightgrey",
             lwd= .75)
  abline(h= 2, lty= "13")
}]
dev.off()