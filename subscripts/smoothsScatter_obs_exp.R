setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="control") | (group_L=="control" & group_R=="dev")), plot_group:= "vllib015: Ctl. x dev / dev x Ctl."]
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="DHS_peak") | (group_L=="DHS_peak" & group_R=="dev")), plot_group:= "vllib015: DHS x dev / dev x DHS"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="dev", plot_group:= "vllib015: dev x dev"]
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="Silencer") | (group_L=="Silencer" & group_R=="dev")), plot_group:= "vllib015: Sil. x dev / dev x Sil."]
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="SUHW_peak") | (group_L=="SUHW_peak" & group_R=="dev")), plot_group:= "vllib015: Put. Ins. x dev / dev x Put. Ins."]
dat[cdition=="vllib016" & group_L=="hk" & group_R=="hk", plot_group:= "vllib016: hk x hk"]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="control") | (group_L=="control" & group_R=="hk")), plot_group:= "vllib016: Ctl. x hk / hk x Ctl."]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="DHS_peak") | (group_L=="DHS_peak" & group_R=="hk")), plot_group:= "vllib016: DHS x hk / hk x DHS"]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="Silencer") | (group_L=="Silencer" & group_R=="hk")), plot_group:= "vllib016: Sil. x hk / hk x Sil."]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="SUHW_peak") | (group_L=="SUHW_peak" & group_R=="hk")), plot_group:= "vllib016: Put. Ins. x hk / dev x hk"]
dat[, plot_group:= factor(plot_group, levels = c("vllib015: Sil. x dev / dev x Sil.",
                                                 "vllib015: Put. Ins. x dev / dev x Put. Ins.",
                                                 "vllib015: Ctl. x dev / dev x Ctl.",
                                                 "vllib015: DHS x dev / dev x DHS",
                                                 "vllib015: dev x dev",
                                                 "vllib016: Sil. x hk / hk x Sil.",
                                                 "vllib016: Put. Ins. x hk / dev x hk", 
                                                 "vllib016: hk x hk", 
                                                 "vllib016: Ctl. x hk / hk x Ctl.", 
                                                 "vllib016: DHS x hk / hk x DHS"))]

# PLOT
pdf("pdf/smoothScatter_observed_expected.pdf", width = 21, height = 9)
par(mfrow= c(2,5))
dat[!is.na(plot_group), {
  lim <- c(-3, switch(cdition, "vllib015"= 12, "vllib016"= 9))
  smoothScatter(additive, 
                log2FoldChange, 
                xlab= "Expected additive (log2)",
                ylab= "Activity (log2)",
                xlim= lim,
                ylim= lim,
                las= 1)
  .lm <- lm(log2FoldChange~additive)
  .eq <- vl_model_equation(.lm, digits= 2)
  mtext(gsub("log2FoldChange", "Activity", .eq), line = 0.4)
  mtext(plot_group, line = 2)
  abline(.lm, lty= 2)
  abline(0, 1)
}, keyby= .(plot_group, cdition)]
dev.off()