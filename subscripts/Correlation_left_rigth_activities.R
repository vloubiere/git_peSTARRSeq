setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds")
dat <- merge(unique(dat[, .(ID= L, indL)]),
             unique(dat[, .(ID= R, indR)]))

# Add TWIST data ----
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
twist <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/BA_300bp_TWIST_STARRSeq.txt")
lib[twist, dev_log2FC_TWIST:= i.dev_log2FoldChange, on= "BA_ID==ID"]
dat[lib, dev_log2FC_TWIST:= i.dev_log2FC_TWIST, on= "ID==ID_vl"]
breaks <- seq(0, nrow(dat), length.out= 5)
col <- c("grey20", "#4D9221", "#B8E186", "#F1B6DA", "#C51B7D")
Cc <- circlize::colorRamp2(breaks, col)
dat[, TWIST_col:= Cc(rank(dev_log2FC_TWIST, na.last = F))]

# Model ----
.lm <- lm(indL~indR, dat)

# Plot ----
pdf("pdf/draft/Correlation_left_rigth_activities.pdf", 
    width = 3, 
    height = 3)
par(mai= c(0.75,0.75,0.75,0.75), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
dat[, {
  # Scatter plot
  plot(indR,
       indL,
       pch= 16,
       cex= 0.5,
       col= adjustcolor(TWIST_col, 0.7),
       las= 1,
       xaxt= "n",
       xlab= "3' activity (log2)",
       ylab= "5' activity (log2)")
  .SD[grepl("^control", ID), {
    points(indR, indL, cex= 0.5, lwd= .25)
  }]
  axis(1, padj = -1.25)
  vl_plot_coeff(value = cor.test(indR, indL)$estimate,
                type= "pcc",
                cex= 7/12)
}]
abline(.lm, lty= "11")
# Legend
vl_heatkey(breaks,
           col,
           left = par("usr")[2]-strwidth("M")*2,
           top = par("usr")[4]-(par("usr")[4]-par("usr")[3])/1.6,
           show.breaks = F,
           main = "Act. rank",
           main.cex = 6/12,
           border= NA,
           width = strwidth("M")/1.5,
           height = strheight("M")*4)
par(lwd= .25)
legend(par("usr")[2]-strwidth("M")*2.3,
       par("usr")[4]-(par("usr")[4]-par("usr")[3])/2.5,
       legend= "Control",
       bty= "n",
       pch= 1,
       cex= 6/12,
       xpd= T)
# # Legend
# unique(dat[, .(TWIST_class, TWIST_col)])[,{
#   legend("bottomright",
#          col= adjustcolor(TWIST_col, 0.7),
#          pch= 16,
#          legend= TWIST_class,
#          bty= "n",
#          cex= 6/12,
#          y.intersp = .75,
#          inset= c(-.3, 0),
#          xpd= T)
# }]
dev.off()