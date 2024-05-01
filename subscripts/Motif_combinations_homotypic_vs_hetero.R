setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

# Import mean act/residuals per motif combination ----
dat <- readRDS("db/motif_effect/motif_impact_act_res.rds")
dat[, V1:= gsub("__L$", "", V1)]
dat[, V2:= gsub("__R$", "", V2)]

# Add names ----
names <- readRDS("db/motif_counts/lib8_motifs_enrichments.rds")
cols <- c("V1", "V2")
dat[, (cols):= lapply(.SD, function(x) gsub("__L$|__R$", "", x)), .SDcols= cols]
dat[names, V1:= i.name, on= "V1==motif_ID"]
dat[names, V2:= i.name, on= "V2==motif_ID"]

# Add colors for pairs of interest ----
dat[V1 %in% c("AP-1", "GATA", "SREBP", "Twist", "Trl") & V1==V2, label:= paste0(V1, "/", V2)]
dat[V1=="AP-1" & V2 %in% c("AP-1", "GATA", "SREBP", "Twist", "Trl"), label:= paste0(V1, "/", V2)]
dat[V1=="GATA" & V2=="SREBP", label:= paste0(V1, "/", V2)]
dat[!is.na(label), col:= ifelse(V1==V2, "tomato", "limegreen")]
dat[, yadj:= 0.01]
dat[, x:= ifelse(V1==V2, 1, 2)]

# Plot ----
pdf("pdf/draft/motif_homotypic_vs_hetero_act_vs_res.pdf", 
    width = 3,
    height = 3)
par(mai = c(.7, 1, .7, .5), 
    mgp= c(1.25, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2)
# Boxplot homotypic vs heterotypic
vl_boxplot(dat[V1==V2, `Mean residuals [motif / no motif] (log2)`],
           dat[V1!=V2, `Mean residuals [motif / no motif] (log2)`],
           compute.pval = list(c(1,2)),
           names= c("Homotypic pairs",
                    "Heterotypic pairs"),
           tilt.names = T,
           ylab= "Mean residuals\n[motif / no motif] (log2)",
           lwd= .75,
           boxwex= .15)
# Add examples
dat[!is.na(col), {
  points(x,
         `Mean residuals [motif / no motif] (log2)`,
         bg= col,
         pch= 21,
         cex= .5,
         lwd= .5)
  segments(x, 
           `Mean residuals [motif / no motif] (log2)`,
           x+.1,
           `Mean residuals [motif / no motif] (log2)`+yadj,
           lwd= .5)
  text(x+.1,
       `Mean residuals [motif / no motif] (log2)`+yadj,
       label,
       pos= 4,
       cex= 5/12,
       offset= .05,
       xpd= T)
}]
dev.off()