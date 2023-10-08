setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- fread("db/peaks/STARR_DSCP_200_peaks.txt")
dat <- dat[`-log10(padj)`>5]
dist <- vl_closestBed(dat, min_dist = 1)

pdf("pdf/draft/enhancer_enhancer_distance.pdf",
    height = 3.75,
    width = 2.25)
par(tcl= -0.2,
    mgp= c(2,.5,0),
    las=1)
vl_boxplot(abs(dist$dist),
           notch= T,
           ylab = "Enh.-enh. distance (kb)",
           yaxt= "n",
           col= "lightgrey")
axis(2, seq(0,14000,2000), seq(0,14000,2000)/1000)
abline(h= 2000, lty= "11")
dev.off()