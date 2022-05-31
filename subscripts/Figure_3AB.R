setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
source("git_peSTARRSeq/functions/plot_transgene.R")

# Import data
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
lib <- readRDS("Rdata/final_results_table.rds")
lib[, diff:= log2FoldChange-additive]
lib <- feat$add_feature(lib, feat$lib)

# Merge dev and housekeeping screens
dat <- merge(lib[vllib=="vllib015" & group_L %in% c("hk", "dev") & group_R %in% c("hk", "dev"), !c("vllib", "CP", "spacer")],
             lib[vllib=="vllib016" & group_L %in% c("hk", "dev") & group_R %in% c("hk", "dev"), !c("vllib", "CP", "spacer")],
             by= c("L", "R", "group_L", "group_R", "class_L", "class_R", "col"), 
             suffixes= c("_dev", "_hk"))
dat[group_L=="dev" & group_R=="dev", c("col", "class"):= .("#74C27A", "dev/dev")]
dat[group_L=="dev" & group_R=="hk", c("col", "class"):= .("cyan", "dev/hk")]
dat[group_L=="hk" & group_R=="dev", c("col", "class"):= .("royalblue2", "hk/dev")]
dat[group_L=="hk" & group_R=="hk", c("col", "class"):= .("tomato", "hk/hk")]
dat[, class:= factor(class, c("dev/dev", "hk/dev", "dev/hk", "hk/hk"))]

pdf("pdf/draft/Figure_3AB.pdf",
    height = 3.25,
    width= 3.3)
layout(matrix(c(1,4,3,2), ncol = 2, byrow = T),
       widths = c(2,0.46),
       heights = c(0.42,2))

for(var in c("log2FoldChange", "diff"))
{
  x <- unlist(dat[, grep(paste0(var, "_dev"), names(dat)), with= F])
  xl <- switch(var, 
               "log2FoldChange"= "Activity dCP (log2)",
               "diff"= "Obs./Exp. Add. dCP (log2)")
  y <- unlist(dat[, grep(paste0(var, "_hk"), names(dat)), with= F])
  yl <- switch(var, 
               "log2FoldChange"= "Activity hkCP (log2)",
               "diff"= "Obs./Exp. Add. hkCP (log2)")
  
  # Activity
  par(mgp= c(1.5, 0.5, 0),
      mar= c(0,4,0,0),
      tcl= -0.2, 
      las= 1,
      xaxs= "r",
      yaxs= "r")
  vl_boxplot(split(x, dat$class),
             horizontal= T,
             compute_pval = list(c(1,4)), 
             violin= T,
             violcol= adjustcolor(dat[, col[1], keyby= class]$V1, 0.5),
             pval_offset = 0.02,
             las= 1,
             xaxt= "n")
  segments(0, 1, 0, 4, 
           lty= 2)
  xlim <- par("usr")[c(1,2)]
  par(mar= c(4,0,0,0))
  vl_boxplot(split(y, dat$class),
             compute_pval = list(c(1,4)), 
             violin= T,
             violcol= adjustcolor(dat[, col[1], keyby= class]$V1, 0.5),
             pval_offset = 0.02,
             yaxt= "n", 
             tilt.names = T)
  
  segments(1, 0, 4, 0, 
           lty= 2)
  ylim <- par("usr")[c(3,4)]
  par(mar= c(4,4,0,0),
      xaxs= "i",
      yaxs= "i")
  plot(x,
       y, 
       xlab= xl,
       ylab= yl,
       xlim= xlim,
       ylim= ylim,
       pch= 16,
       col= adjustcolor(dat$col, 0.5),
       cex= 0.5,
       xaxt= "n",
       yaxt= "n")
  axis(1,
       axisTicks(xlim, 
                 log= F, 
                 nint = 4),
       lwd= 0,
       lwd.ticks= 1)
  axis(2,
       axisTicks(ylim, 
                 log= F, 
                 nint = 4),
       lwd= 0,
       lwd.ticks= 1)
  abline(h= 0, lty= 2)
  abline(v= 0, lty= 2)
  abline(0,1,lty= 2)
  legend("topleft",
         pch= 16,
         legend= dat[, class[1], keyby= class]$V1,
         col= adjustcolor(dat[, col[1], keyby= class]$V1, 0.8),
         bty= "n",
         cex= 0.7)
  # Transgenes
  transgene(x= par("usr")[1]-diff(grconvertX(c(0,3), "lines", "user")),
            y= grconvertY(0.435, "nfc", "user"), 
            CP = "hk", 
            rotate= T,
            cex= 0.6)
  transgene(x= grconvertX(0.435, "nfc", "user"),
            y= par("usr")[3]-diff(grconvertY(c(0,3.25), "lines", "user")), 
            CP = "dev",
            cex= 0.6)
  par(mar= c(0,0,0,0))
  plot.new()
}
dev.off()