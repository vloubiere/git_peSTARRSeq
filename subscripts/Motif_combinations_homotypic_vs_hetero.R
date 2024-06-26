setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

# Import mean act/residuals per motif combination ----
dat <- readRDS("db/motif_effect/motif_impact_act_res.rds")
dat[, V1:= gsub("__L$", "", V1)]
dat[, V2:= gsub("__R$", "", V2)]

# Select combinations of interest ----
sel <- data.table(name= c("GATA", "AP-1", "Twist", "Trl", "SREBP"),
                  motif_ID= c("cisbp__M4320",
                              "cisbp__M6317",
                              "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715",
                              "homer__CTCTCTCTCY_GAGA-repeat",
                              "cisbp__M2388"))
cmb <- data.table(V1= c("GATA", "AP-1", "Twist", "Trl", "SREBP", "AP-1", "GATA", "AP-1", "GATA"),
                  V2= c("GATA", "AP-1", "Twist", "Trl", "SREBP", "Twist", "Trl", "GATA", "SREBP"))
cmb[, label:= paste0(V1, "/", V2)]
cmb <- merge(cmb, sel, by.x= "V1", by.y= "name", all.x= T)
cmb <- merge(cmb, sel, by.x= "V2", by.y= "name", all.x= T)
cmb$V1 <- cmb$V2 <- NULL

# Merge to data ----
dat <- merge(dat, cmb, by.x= c("V1", "V2"), by.y= c("motif_ID.x", "motif_ID.y"), all.x= T)
dat[!is.na(label), col:= ifelse(V1==V2, "tomato", "limegreen")]
dat[, yadj:= 0.01]
dat[, x:= ifelse(V1==V2, 1, 2)]

# Plot ----
pdf("pdf/draft/motif_homotypic_vs_hetero_act_vs_res.pdf", 
    width = 3.15,
    height = 3)
par(mai = rep(.9, 4), 
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
},.(x, `Mean residuals [motif / no motif] (log2)`, col, yadj)]
dev.off()