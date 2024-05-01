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
dat[names, V1:= i.name, on= "V1==motif_ID"]
dat[names, V2:= i.name, on= "V2==motif_ID"]

# Add colors for pairs of interest ----
dat[, col:= fcase(V1 %in% c("Twist", "Trl") & V2 %in% c("Twist", "Trl"), "tomato",
                  V1 %in% c("Dref.1", "Dref.2") & V2 %in% c("Dref.1", "Dref.2"), "cornflowerblue",
                  V1 %in% c("AP-1", "GATA") & V2 %in% c("AP-1", "GATA"), "limegreen",
                  default= "lightgrey")]
dat[col=="limegreen", col:= ifelse(`Mean activity [motif / no motif] (log2)`<max(`Mean activity [motif / no motif] (log2)`), "lightgrey", col), .(V1, V2)]
dat[col!="lightgrey", label:= paste0(V1, "/", V2)]
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