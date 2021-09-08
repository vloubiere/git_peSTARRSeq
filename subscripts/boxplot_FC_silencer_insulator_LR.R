setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")
sel <- c("vllib015: dev x Silencer / Silencer x dev",
         "vllib015: dev x SUHW_peak / SUHW_peak x dev",
         "vllib016: hk x Silencer / Silencer x hk",
         "vllib016: hk x SUHW_peak / SUHW_peak x hk")
dat <- DT[plot_group %in% sel & L!=R]
dat[, plot_group:= factor(plot_group, levels = sel)]
dat <- dat[!is.na(plot_group)]
dat[group_L %in% c("dev", "hk"), FC:= log2FoldChange-median_L]
dat[group_R %in% c("dev", "hk"), FC:= log2FoldChange-median_R]

# PLOT
pdf("pdf/aggregate_activity/boxplot_FC_sil_insulator_LR.pdf", 
    width = 10,
    height= 4.5)
par(las= 1, 
    mfrow= c(1,2))
dat[, {
  .c <- split(FC, plot_group_LR)
  names(.c) <- gsub(paste0(cdition, ": "), "", names(.c))
  vl_boxplot(.c, 
             col= "black", 
             ylim = c(-6, 4))
  abline(h= 0, 
         lty= 2)
  title(cdition)
  segments(1,3.2,2,3.2)
  vl_plot_pval(x = 1.5,
               y= 3.6,
               wilcox.test(.c[[1]], .c[[2]])$p.value, 
               cex= 0.8)
  print("")
}, .(plot_group, cdition)]
dev.off()

