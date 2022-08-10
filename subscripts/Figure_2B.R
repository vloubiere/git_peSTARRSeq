setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import Clustering object (see clustering_vllib002_lm_residuals_Figure_2B.R)
cl <- readRDS("Rdata/vllib002_clustering_expected_scores_draft_figure.rds")
ind_L <- cl$rows$median
ind_R <- cl$cols$median

# plot
pdf("pdf/draft/Figure_2B.pdf", 
    width= 8, 
    height = 7)
par(mar= c(3.5,3.5,3,10),
    mgp= c(3,0.15,0),
    tcl= -0.2)

# Heatmap
plot(cl)

# Left individual activities
right <- par("usr")[1]-strwidth("M", cex= 0.5)
width <- right-grconvertX(1, "line", "user")
rect(right-(ind_L/max(ind_L)*width), 
     rev(seq(ind_L))+0.5, 
     right, 
     rev(seq(ind_L))-0.5, 
     border= NA, 
     xpd= T, 
     col= adjustcolor("#0C3A0E", 0.7))
ticks <- axisTicks(c(0, max(ind_L)), log= F)
at <- right-(ticks/max(ind_L)*width)
axis(3, 
     at = rev(range(at)), 
     labels = range(ticks),
     xpd= T, 
     line= 0.25, 
     cex.axis= 0.5)
text(mean(at),
     par("usr")[4],
     "5' Individual\nact. (log2)",
     xpd= T, 
     pos= 3, 
     offset= 1.25,
     cex= 0.7)
text(grconvertX(0.5, "line", "user"),
     mean(par("usr")[c(3,4)]),
     "5' enhancer",
     srt= 90,
     xpd= T)

# Right individual activities
par(mgp= c(3,0.35,0))
bottom <- grconvertY(1, "line", "user")
height <- par("usr")[3]-strheight("M")*0.5-bottom
rect(seq(ind_R)+0.5,
     bottom+(ind_R/max(ind_R)*height), 
     seq(ind_R)-0.5,
     bottom,
     border= NA, 
     xpd= T, 
     col= adjustcolor("#0C3A0E", 0.7))
ticks <- axisTicks(c(0, max(ind_R)), log= F)
at <- bottom+(ticks/max(ind_R)*height)
axis(2, 
     at = range(at), 
     labels = range(ticks),
     xpd= T, 
     line= 0.25, 
     cex.axis= 0.5,
     las= 2)
text(par("usr")[1],
     mean(at),
     "3' Individual\nact. (log2)",
     xpd= T,
     pos= 2, 
     offset= 0.5,
     cex= 0.7)
text(mean(par("usr")[c(1,2)]),
     grconvertY(0.5, "line", "user"),
     "3' enhancer",
     xpd= T)
dev.off()