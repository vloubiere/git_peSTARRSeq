setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="control") | (group_L=="control" & group_R=="dev")), plot_group:= "vllib015: Ctl. x dev / dev x Ctl."]
dat[cdition=="vllib015" & ((group_L=="dev" & group_R=="SUHW_peak") | (group_L=="SUHW_peak" & group_R=="dev")), plot_group:= "vllib015: Put. Ins. x dev / dev x Put. Ins."]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="control") | (group_L=="control" & group_R=="hk")), plot_group:= "vllib016: Ctl. x hk / hk x Ctl."]
dat[cdition=="vllib016" & ((group_L=="hk" & group_R=="SUHW_peak") | (group_L=="SUHW_peak" & group_R=="hk")), plot_group:= "vllib016: Put. Ins. x hk / hk x Put. Ins."]

# PLOT
pdf("pdf/boxplot_insulator_act_LR.pdf", width = 15, height = 5)
dat[!is.na(plot_group), {
  left_enh_ins <- grepl("Put. Ins.", plot_group) & group_L==ifelse(cdition=="vllib015", "dev", "hk")
  right_enh_ins <- grepl("Put. Ins.", plot_group) & group_R==ifelse(cdition=="vllib015", "dev", "hk")
  left_enh_ctl <- grepl("Ctl.", plot_group) & group_L==ifelse(cdition=="vllib015", "dev", "hk")
  right_enh_ctl <- grepl("Ctl.", plot_group) & group_R==ifelse(cdition=="vllib015", "dev", "hk")
  # Compute left/right and ctl/ins pval
  pval <- c(wilcox.test(log2FoldChange[left_enh_ins], 
                        log2FoldChange[right_enh_ins])$p.value,
            wilcox.test(log2FoldChange[left_enh_ctl],
                        log2FoldChange[right_enh_ctl])$p.value,
            wilcox.test(log2FoldChange[left_enh_ins],
                        log2FoldChange[left_enh_ctl])$p.value,
            wilcox.test(log2FoldChange[right_enh_ins],
                        log2FoldChange[right_enh_ctl])$p.value)
  # Compute left/right and ctl/ins FC
  diff <- c(median(log2FoldChange[right_enh_ins])/median(log2FoldChange[left_enh_ins]),
            median(log2FoldChange[right_enh_ctl])/median(log2FoldChange[left_enh_ctl]),
            median(log2FoldChange[left_enh_ctl])/median(log2FoldChange[left_enh_ins]),
            median(log2FoldChange[right_enh_ctl])/median(log2FoldChange[right_enh_ins]))
  diff <- paste0("x", round(diff, 2))
  # Boxplot
  boxplot(log2FoldChange[left_enh_ins],
          log2FoldChange[right_enh_ins],
          log2FoldChange[left_enh_ctl],
          log2FoldChange[right_enh_ctl],
          names= NA,
          las= 1,
          ylab= "Activity (log2)",
          whisklty= 1,
          staplelwd = NA,
          ylim= c(-4, 12.5),
          pch= 19,
          cex= 0.4,
          main= cdition)
  # Boxplot labels
  labs <- c("DEV. / Su(H)w", "Su(H)w / DEV.", "DEV. / Control", "Control / DEV.")
  if(cdition=="vllib016")
    labs <- gsub("DEV", "HK", labs)
  text(1:4+strwidth("D")*1.5,
       grconvertY(0, "npc", "user")-strheight("D")*1.5,
       labs,
       srt= 45,
       pos= 2,
       xpd= T)
  # Plot pval
  segments(c(1,3,1,2),
           c(10, 10, 11, 12),
           c(2,4,3,4),
           c(10, 10, 11, 12))
  vl_plot_pval(x = c(1.5, 3.5, 2, 3)+strwidth(diff)/2,
               y = c(10,10,11,12)+grconvertY(strheight("a"))*0.6,
               pval = pval,
               stars_only = T,
               pos= 4,
               offset= 0)
  # Plot diff
  text(x = c(1.5, 3.5, 2, 3),
       y = c(10, 10, 11, 12)+grconvertY(strheight("a"))*0.6,
       diff,
       cex= 0.9)
  abline(h= max(c(median(log2FoldChange[left_enh_ins]),
                  median(log2FoldChange[right_enh_ins]))), 
         lty= 2)
  print("")
}, cdition]
dev.off()


