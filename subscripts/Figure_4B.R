setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table_CP_basalAct.rds")
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
dat[, CP:= factor(CP, c("dev", "hk"))]
dat[, basalAct:= factor(basalAct, c("ref", "low", "high"))]

pdf("pdf/draft/Figure_4B.pdf",
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
      if(CP=="hk")
        margins[c(1, 3)] <- c(6, 0.05)
      par(mar= margins,
          mgp= c(2, 0.5, 0),
          lwd= 0.5)
      vl_boxplot(diff~basalAct,
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
      if(CP=="dev")
        title(class)
      if(CP=="hk")
        axis(1, 
             at= seq(unique(basalAct)),
             labels = unique(sort(basalAct)),
             las= 2,
             tck= -0.025,
             lwd= 0.5)
      if(class=="hk/hk")
        title(ylab= "Observed/Exp. add.", line = )
      print("")
    }, keyby= .(CP, class, Cc)]
dev.off()
