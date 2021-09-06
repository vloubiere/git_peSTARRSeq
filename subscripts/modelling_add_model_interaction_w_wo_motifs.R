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

#------------------------#
# Build models 
#------------------------#
colnames(dat) <- gsub("-", "_", colnames(dat))

pdf("pdf/modeling/additive_model_with_interaction.pdf", width = 13)
par(mfrow= c(1,2))
dat[!is.na(plot_group) & !is.na(median_L) & !is.na(median_R) & L!=R, {
  # Additive model
  .lm1 <- lm(log2FoldChange~median_L*median_R)
  pred1 <- predict(.lm1)
  pcc1 <- round(cor.test(log2FoldChange, pred1)$estimate, 2)
  smoothScatter(pred1, 
                log2FoldChange,
                ylab= "Activity (log2)",
                xlab= "Predicted (log2)",
                las= 1)
  mtext(plot_group)
  abline(0, 1, lty= 2)
  legend("topleft", 
         legend = paste0("PCC= ", pcc1), 
         bty= "n")
  # With motifs
  .form <- paste0("log2FoldChange~median_L*median_R+", paste0(grep("motif", names(.SD), value = T), collapse = "+"))
  .lm2 <- lm(as.formula(.form))
  pred2 <- predict(.lm2)
  pcc2 <- round(cor.test(log2FoldChange, pred2)$estimate, 2)
  smoothScatter(pred2,
                log2FoldChange,
                ylab= "Activity (log2)",
                xlab= "Predicted (log2)",
                las= 1)
  mtext(plot_group)
  abline(0, 1, lty= 2)
  legend("topleft",
         legend = paste0("PCC= ", pcc2),
         bty= "n")
  print("")
}, plot_group]
dev.off()
