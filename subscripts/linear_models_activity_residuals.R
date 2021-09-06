setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
if("plot_group" %in% names(dat))
  dat$plot_group <- NULL
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="control" & median_L>1) | (group_L=="control" & group_R=="dev" & median_R>1)), plot_group:= "vllib002: Ctl. x dev / dev x Ctl."]
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="hk" & median_L>1) | (group_L=="hk" & group_R=="dev" & median_R>1)), plot_group:= "vllib002: dev x hk / dev x hk"]
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x dev"]
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
dat[, plot_group:= factor(plot_group, levels = c("vllib002: Ctl. x dev / dev x Ctl.",
                                                 "vllib002: inducible x dev / dev x inducible",
                                                 "vllib002: OSC x dev / dev x OSC",
                                                 "vllib002: dev x hk / dev x hk",
                                                 "vllib002: dev x dev",
                                                 "vllib015: Sil. x dev / dev x Sil.",
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
pdf("pdf/aggregate_activity/smoothScatter_observed_expected.pdf", 
    width = 18*0.7, 
    height = 4.7*0.7)
layout(matrix(1:5, ncol= 5), widths = c(1, rep(0.82, 3), 1))
dat[!is.na(plot_group) & L!=R, {
  lim <- c(switch(cdition, 
                  "vllib002"= -1, 
                  "vllib015"= -3, 
                  "vllib016"= -3),
           switch(cdition, 
                  "vllib002"= 10, 
                  "vllib015"= 12, 
                  "vllib016"= 9))
  if(.GRP %in% seq(1,101, 5))
  {
    pl <- T
    yl <- "Activity (log2)"
    par(mar= c(5.1, 4.1, 4.1, 0))
  }else if(.GRP %in% seq(5,100, 5))
  {
    pl <- F
    yl <- NA
    par(mar= c(5.1, 0, 4.1, 4.1))
  }else
  {
    pl <- F
    yl <- NA
    par(mar= c(5.1, 0, 4.1, 0))
  }
  smoothScatter(additive, 
                log2FoldChange, 
                xlab= NA,
                ylab= NA,
                xlim= lim,
                ylim= lim,
                yaxt= ifelse(pl, "s", "n"),
                las= 1)
  .lm <- lm(log2FoldChange~additive)
  .eq <- vl_model_equation(.lm, digits= 2)
  mtext(side=1, 
        text= "Expected additive (log2)", 
        line= 2.5,
        cex= 0.8)
  mtext(side=2.5, 
        text= yl, 
        line= 2,
        cex= 0.8)
  mtext(gsub("log2FoldChange", "Activity", .eq), 
        line = 0.8,
        cex= 0.8)
  mtext(plot_group, 
        line = 3,
        cex= 0.8)
  abline(.lm, lty= 2)
  abline(0, 1)
}, keyby= .(plot_group, cdition)]
dev.off()
