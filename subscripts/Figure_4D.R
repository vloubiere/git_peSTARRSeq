setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(data.table)
source("git_peSTARRSeq/functions/plot_transgene.R")

#---------------------------------------------#
# Import full data
#---------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")

#---------------------------------------------#
# Format
#---------------------------------------------#
dat <- lib[vllib %in% c("vllib015", "vllib017", "vllib023", "vllib018")
           & class=="enh./enh." 
           & class=="enh./enh." 
           & grepl("^dev", L) 
           & grepl("^dev", R)][, diff:= log2FoldChange-additive]
dat <- dat[, check:= all(unique(dat$spacer) %in% spacer), .(L, R)][(check)]
dat[, spacer:= switch(spacer,
                      "SCR1_300"= "no intron",
                      "intron2_300"= "300bp intron1",
                      "shortened-intron4_300"= "300bp intron2",
                      "intron4_2000"= "2kb intron2"), spacer]
dat[, spacer:= factor(spacer, 
                      c("no intron",  "300bp intron1", "300bp intron2", "2kb intron2"))]
ind <- rbind(dat[, .(enh= L, act= median_L, spacer, side= "5'")],
             dat[, .(enh= R, act= median_R, spacer, side= "3'")])
ind <- unique(ind)
ind[, side:= factor(side, c("5'", "3'"))]

pdf("pdf/draft/Figure_4D.pdf",
    height= 3,
    width = 4)
layout(matrix(1:2, ncol= 2),
       widths = c(1, 0.4))
par(mar= c(3,3,1,4),
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2)
vl_boxplot(act~spacer+side, 
           ind, 
           compute_pval= list(c(1,2), c(3,4)),
           ylab= "Individual activity (log2)",
           notch= T,
           boxcol= adjustcolor(c("#33FF99", "#FFCCFF"), 0.6),
           las= 2)
segments(1, 
         par("usr")[3], 
         2, 
         par("usr")[3], xpd= T)
text(1.5, 
     par("usr")[3], 
     "5'", 
     pos= 1,
     xpd= T)
segments(3, 
         par("usr")[3], 
         4, 
         par("usr")[3], xpd= T)
text(3.5, 
     par("usr")[3], 
     "3'", 
     pos= 1,
     xpd= T)
legend(4,
       par("usr")[4],
       legend = c("Low basal act.",
                  "High basal act."),
       fill= adjustcolor(c("#33FF99", "#FFCCFF"), 0.6),
       bty= "n",
       xpd= T,
       cex= 0.7)
par(mar= c(3,2.5,1,0))
vl_boxplot(diff~spacer,
           dat,
           compute_pval= list(c(1,2)), 
           tilt.names= T,
           ylab= "Obs./Exp. Add. (log2)",
           notch= T,
           boxcol= adjustcolor(c("#33FF99", "#FFCCFF"), 0.6),
           ylim= c(-1,4.5))
abline(h=0, 
       lty=2)
dev.off()
