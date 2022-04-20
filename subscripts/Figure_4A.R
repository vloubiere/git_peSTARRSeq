setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table_spacer_size.rds")
dat <- dat[act_wilcox_L<0.001 & median_L>log2(1.5) 
           & act_wilcox_R<0.001 
           & median_R>log2(1.5) 
           & grepl("^dev|^hk", L) 
           & grepl("^dev|^hk", R)]
dat[, diff:= log2FoldChange-additive]
set.seed(1)
dat <- dat[sample(nrow(dat), nrow(dat))]
dat[, class:= factor(
  fcase(grepl("dev", L) & grepl("dev", R), "dev/dev",
        grepl("dev", L) & grepl("hk", R), "dev/hk",
        grepl("hk", L) & grepl("dev", R), "hk/dev",
        grepl("hk", L) & grepl("hk", R), "hk/hk"),
  levels= c("hk/hk",
            "hk/dev",
            "dev/hk",
            "dev/dev"))]
dat[class=="dev/dev", Cc:= "#74C27A"]
dat[class=="hk/dev", Cc:= "royalblue2"]
dat[class=="dev/hk", Cc:= "cyan"]
dat[class=="hk/hk", Cc:= "tomato"]
dat[, CP:= factor(CP, c("DSCP", "RpS12"))]
dat[, spacer:= factor(spacer, c("no_intron", "short_intron", "long_intron"))]

pdf("pdf/draft/Figure_4A.pdf",
    height= 3,
    width = 4)
layout(matrix(1:8, 
              nrow= 2, 
              byrow = T),
       widths = c(1.3,1,1,1),
       heights = c(1, 1.375))
dat[,
    {
      margins <- c(1,2,1.5,0.1)
      if(class=="hk/hk")
        margins[2] <- 4
      if(CP=="RpS12")
        margins[c(1, 3)] <- c(6, 0.05)
      par(mar= margins,
          mgp= c(2, 0.5, 0),
          lwd= 0.5)
      vl_boxplot(diff~spacer,
                 ylim= c(-4.5, 5.5),
                 boxcol= Cc,
                 xaxt= "n",
                 yaxt= "n",
                 compute_pval= list(c(1,2), c(2,3)),
                 notch= T)
      axis(2, 
           at= seq(-4, 4, 2),
           labels = seq(-4, 4, 2),
           las= 2,
           tck= -0.025,
           lwd= 0.5)
      abline(h= 0,
             lty= 2)
      if(CP=="DSCP")
        title(class)
      if(CP=="RpS12")
        axis(1, 
             at= seq(unique(spacer)),
             labels = unique(sort(spacer)),
             las= 2,
             tck= -0.025,
             lwd= 0.5)
      if(class=="hk/hk")
        title(ylab= "Observed/Exp. add.", line = )
      print("")
    }, keyby= .(CP, class, Cc)]
dev.off()

# pdf("pdf/draft/Figure_4A.pdf", 
#     width = 7.5)
# mat <- matrix(1:16, 
#               ncol= 4, 
#               byrow = T)
# layout(mat, 
#        heights = c(1,1,1,0.25),
#        widths = c(1,0.5,1,0.5))
# dat[, {
#   # Scatterplot
#   par(mar= c(2,4,0.1,0.1),
#       las= 1,
#       mgp = c(1.5, 0.5, 0))
#   if(CP=="DSCP")
#     lim <- c(-1.75, 11) else if(CP=="RpS12")
#       lim <- c(-3, 10)
#   plot(NA, 
#        xlim= lim,
#        ylim= lim,
#        ylab= "Activity (log2)",
#        xlab= NA,
#        xaxt= "n",
#        yaxt= "n")
#   axis(1,
#        at= seq(-2, 10, 2),
#        labels= seq(-2, 10, 2),
#        tck= -0.025)
#   axis(2,
#        at= seq(-2, 10, 2),
#        labels= seq(-2, 10, 2),
#        tck= -0.025)
#   points(additive, 
#          log2FoldChange,
#          col= adjustcolor(Cc, 0.5),
#          pch= 19,
#          cex= 0.5)
#   legend("topleft", 
#          legend= "Super-additive",
#          bty= "n")
#   legend("bottomright", 
#          legend= "Sub-additive",
#          bty= "n")
#   abline(0, 1)
#   # Axes
#   if(.GRP %in% c(5, 6))
#     title(xlab= "Expected additive (log2)", 
#           xpd= T)
#   # Boxplot per class
#   par(mar= c(2,3,0.1,2))
#   vl_boxplot(diff~class, 
#              yaxt= "n",
#              xaxt= "n", 
#              boxcol= unique(Cc[order(class)]))
#   axis(2, 
#        axisTicks(range(diff), log = F),
#        axisTicks(range(diff), log = F), 
#        tck= -0.075)
#   title(ylab= "Observed/Exp. add.")
#   abline(h=0, 
#          lty= 2)
#   # Axes
#   if(.GRP %in% c(5, 6))
#   {
#     axis(1, 
#          seq(unique(class)),
#          unique(sort(class)), 
#          tck= -0.075,
#          las= 2)
#   }
#   print("DONE")
# }, keyby= .(spacer, CP)]
# dev.off()
