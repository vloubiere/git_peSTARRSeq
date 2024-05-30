setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/linear_models/FC_DSCP_long_spacer_act_pairs_lm_predictions.rds")
short <- readRDS("db/FC_tables/DSCP_large_WT_FC_DESeq2.rds")
dat[short, indLshort:= i.indL, on= "L"]
dat[short, indRshort:= i.indR, on= "R"]

# Individual act ----
indL <- unique(dat[, .(L, indL, indLshort)])
indR <- unique(dat[, .(R, indR, indRshort)])

# Plot ----
Cc <- c("white", "lightgrey")

pdf("pdf/draft/long_spacer_individual_act.pdf", 2.25, 3)
vl_par(mai= c(.9, .6, .9, .9),
       mgp= c(.9, .35, 0),
       font.main= 1)
# Boxplot individual act
vl_boxplot(indL$indLshort,
           indL$indL,
           indR$indRshort,
           indR$indR,
           compute.pval = list(c(1,2), c(3,4)),
           ylab= "Activity (log2)",
           xaxt= "n",
           at= c(0.9, 1.1, 1.9, 2.1),
           boxwex= .15,
           col= Cc,
           lwd= .75,
           main= "2kb spacer")
abline(h= 0, lty= "13")
vl_tilt_xaxis(1:2,
              labels = c("5' enhancer", "3' enhancer"))
legend(par("usr")[2]-strwidth("M"),
       par("usr")[4],
       fill= Cc,
       legend = c("300bp spacer", "2kb spacer"),
       bty= "n",
       cex= .7,
       xpd= T)
dev.off()