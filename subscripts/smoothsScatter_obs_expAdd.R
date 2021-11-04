setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")
sel <- c("vllib002: control x dev / dev x control",
         "vllib002: dev x inducible / inducible x dev",
         "vllib002: dev x OSC / OSC x dev",
         "vllib002: dev x hk / hk x dev",
         "vllib002: dev x dev",
         "vllib015: dev x Silencer / Silencer x dev",
         "vllib015: dev x SUHW_peak / SUHW_peak x dev",
         "vllib015: control x dev / dev x control",
         "vllib015: dev x DHS_peak / DHS_peak x dev",
         "vllib015: dev x dev",
         "vllib016: hk x Silencer / Silencer x hk",
         "vllib016: hk x SUHW_peak / SUHW_peak x hk", 
         "vllib016: hk x hk", 
         "vllib016: control x hk / hk x control",
         "vllib016: DHS_peak x hk / hk x DHS_peak")
dat <- DT[plot_group %in% sel & L!=R]
dat <- dat[cdition!="vllib002" | (median_L>1 & median_R>1)]
dat[, plot_group:= factor(plot_group, levels = sel)]
dat <- dat[!is.na(plot_group)]

# PLOT
pdf("pdf/aggregate_activity/smoothScatter_observed_expected.pdf", 
    width = 18*0.7, 
    height = 4.7*0.7)
layout(matrix(1:5, ncol= 5), widths = c(1, rep(0.82, 3), 1))
dat[, {
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
