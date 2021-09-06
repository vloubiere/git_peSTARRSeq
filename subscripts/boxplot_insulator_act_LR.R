setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")
sel <- c("vllib015: control x dev / dev x control",
         "vllib015: dev x SUHW_peak / SUHW_peak x dev",
         "vllib016: control x dev / dev x control",
         "vllib016: dev x SUHW_peak / SUHW_peak x dev")
dat <- DT[plot_group %in% sel & L!=R]
dat[, plot_group:= factor(plot_group, levels = sel)]
dat <- dat[!is.na(plot_group)]


# PLOT
pdf("pdf/aggregate_activity/boxplot_insulator_act_LR.pdf", width = 15, height = 5)
dat[, {
  left_enh_ins <- grepl("SUHW_peak", plot_group) & group_L==ifelse(cdition=="vllib015", "dev", "hk")
  right_enh_ins <- grepl("SUHW_peak", plot_group) & group_R==ifelse(cdition=="vllib015", "dev", "hk")
  left_enh_ctl <- grepl("control", plot_group) & group_L==ifelse(cdition=="vllib015", "dev", "hk")
  right_enh_ctl <- grepl("control", plot_group) & group_R==ifelse(cdition=="vllib015", "dev", "hk")
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
}, plot_group]
dev.off()


