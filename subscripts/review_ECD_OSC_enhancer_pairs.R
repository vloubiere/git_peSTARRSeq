setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- list("S2"= readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds"),
            "S2+ECD"= readRDS("db/FC_tables/DSCP_ECD_WT_DESeq2.rds"),
            "OSC"= readRDS("db/FC_tables/DSCP_OSC_WT_DESeq2.rds"))
dat <- rbindlist(dat, idcol = "cdition")
dat[, cdition:= factor(cdition, c("S2", "S2+ECD", "OSC"))]
dat[, classL:= tstrsplit(L, "_", keep= 1), L]
dat[, classR:= tstrsplit(R, "_", keep= 1), R]

# Enhancer pairs ----
pl <- dat[classL %in% c("control", "ecdysone", "OSC")]
pl <- pl[classL==classR]
pl[, name:= fcase(classL=="control", "Ctl./Ctl.",
                   classL=="ecdysone", "Ecd./Ecd.",
                   classL=="OSC", "OSC/OSC")]
pl[, name:= factor(name, c("Ctl./Ctl.", "Ecd./Ecd.", "OSC/OSC"))]

# Plot ----
Cc <- c("white", "cadetblue2", "tan1")
x <- seq(levels(pl$name))

pdf("pdf/review_ecd_osc_enh_pairs_act.pdf", 3.75, 3)
vl_par(lwd= .5)
vl_boxplot(log2FoldChange~cdition+name,
           pl,
           at= rep(x, each= 3)+c(-0.15, 0, 0.15),
           boxwex= .125,
           xaxt= "n",
           compute.pval= list(c(1,2), c(1,3), c(4,5), c(7,9)),
           col= Cc,
           ylab= "Activity (log2)")
legend(par("usr")[2],
       par("usr")[4],
       xpd= T,
       fill= Cc,
       legend= levels(pl$cdition),
       bty= "n",
       cex= 7/12)
abline(h= 0, lty= 2)
axis(1,
     at= x,
     labels = levels(pl$name),
     lwd= 0)
dev.off()