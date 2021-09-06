setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
# dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="control") | (group_L=="control" & group_R=="dev")), plot_group:= "vllib015: Ctl. x dev / dev x Ctl."]
# dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="DHS_peak") | (group_L=="DHS_peak" & group_R=="dev")), plot_group:= "vllib015: DHS x dev / dev x DHS"]
# dat[cdition=="vllib015" & group_L=="dev" & group_R=="dev", plot_group:= "vllib015: dev x dev"]
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="Silencer") | (group_L=="Silencer" & group_R=="dev")), plot_group:= "vllib015: Sil. x dev / dev x Sil."]
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="SUHW_peak") | (group_L=="SUHW_peak" & group_R=="dev")), plot_group:= "vllib015: Put. Ins. x dev / dev x Put. Ins."]
# dat[cdition=="vllib016" & group_L=="hk" & group_R=="hk", plot_group:= "vllib016: hk x hk"]
# dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="control") | (group_L=="control" & group_R=="hk")), plot_group:= "vllib016: Ctl. x hk / hk x Ctl."]
# dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="DHS_peak") | (group_L=="DHS_peak" & group_R=="hk")), plot_group:= "vllib016: DHS x hk / hk x DHS"]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="Silencer") | (group_L=="Silencer" & group_R=="hk")), plot_group:= "vllib016: Sil. x hk / hk x Sil."]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="SUHW_peak") | (group_L=="SUHW_peak" & group_R=="hk")), plot_group:= "vllib016: Put. Ins. x hk / dev x hk"]
# dat[, plot_group:= factor(plot_group, levels = c("vllib015: Sil. x dev / dev x Sil.",
#                                                  "vllib015: Put. Ins. x dev / dev x Put. Ins.",
#                                                  "vllib015: Ctl. x dev / dev x Ctl.",
#                                                  "vllib015: DHS x dev / dev x DHS",
#                                                  "vllib015: dev x dev",
#                                                  "vllib016: Sil. x hk / hk x Sil.",
#                                                  "vllib016: Put. Ins. x hk / dev x hk", 
#                                                  "vllib016: hk x hk", 
#                                                  "vllib016: Ctl. x hk / hk x Ctl.", 
#                                                  "vllib016: DHS x hk / hk x DHS"))]

# PLOT
pdf("pdf/aggregate_activity/boxplot_FC_sil_insulator_LR.pdf", 
    width = 8*1.5,
    height= 3*1.5)
par(las= 1, 
    mfrow= c(1,2))
dat[!is.na(plot_group) & L!=R, 
{
  if(grepl("Sil.", plot_group))
  {
    class <- "Silencer"
    names <- c("Enhancer\nx\nSilencer", "Silencer\nx\nEnhancer")
  }
  if(grepl("Ins.", plot_group))
  {
    class <- "SUHW_peak"
    names <- c("Enhancer\nx\nInsulator", "Insulator\nx\nEnhancer")
  }
  .c <- list(.SD[group_R==class, log2FoldChange-median_L], 
             .SD[group_L==class, log2FoldChange-median_R])
  vl_boxplot(.c,
             col= "black", 
             plot_labels = F, 
             ylab= paste0("FC ", class, "/Control"), 
             ylim= c(-6, 3))
  abline(h= 0, 
         lty= 2)
  title(cdition)
  text(1:2, 
       grconvertY(0, "npc", "user"),
       pos= 1,
       offset= 1.5,
       names,
       xpd= T)
}, keyby= .(plot_group, cdition)]
dev.off()
