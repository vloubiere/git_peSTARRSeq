setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
options(datatable.print.topn= 1)
require(data.table)
require(pheatmap)
require(circlize)
require(seqLogo)
require(png)
require(colorspace)

feat <- readRDS("Rdata/library/lib_features.rds")
dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[enh_L %in% feat[group %in% c("dev", "shared"), ID] &
             enh_R %in% feat[group %in% c("dev", "shared"), ID] & 
             median_L>1 & median_R>1 & !is.na(diff) & !is.na(median_L) & !is.na(median_R)]
feat <- feat[ID %in% c(dat$enh_L, dat$enh_R)]
setkeyv(feat, "ID")
mot <- readRDS("Rdata/modeling_peSTARRSeq/motifs_linear_models_results.rds")
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")


sel <- unique(mot[, mot1_s:Dmel2])
sel[grepl("^18.", Dmel1) & grepl("^24.fru", Dmel2), effect:= "pos"]
sel[grepl("^30.", Dmel1) & grepl("^26.", Dmel2), effect:= "neg"]
sel <- sel[!is.na(effect)]
sel <- sel[, dat[, .(enh_L, enh_R, diff, log2FoldChange, median_L, median_R)], (sel)]
sel[, count1:= feat[enh_L][[col]]>0, .(col= paste0("motif__", mot1_s))]
sel[, count2:= feat[enh_R][[col]]>0, .(col= paste0("motif__", mot2_s))]
sel <- sel[between(log2FoldChange, 4, 6) & median_L<4 & median_R<4]
sel[, all:= "all"]
sel <- sel[order(factor(effect, levels = c("pos", "neg")))]
ex <- max(sel[, .GRP, mot1_s:mot2_s]$GRP)
  
pdf("pdf/peSTARRSeq/examples_lm_motifs_LR_interaction.pdf", 7*ex, 7)
lay <- matrix(c(1,2,6,3,6,4,6,5), ncol= 2, byrow= T)
for(i in 2:ex)
{
  lay <- cbind(lay, lay+max(lay))
}
layout(lay, widths = rep(c(0.25, 1), ex), heights = rep(c(1.3,1,1,1), ex))
par(las= 1)
sel[, 
    {
      # All
      par(mar= c(3.5,5,4,2))
      my_boxplot(diff~all, .SD)
      mtext("all", side = 1, line = 1)
      
      # Residuals
      Cc <- ifelse(effect=="pos", "tomato", "cornflowerblue")
      box <- my_boxplot(diff~count1+count2, .SD, col_box = Cc)
      text(grconvertX(0.25, "npc"), grconvertY(1.32, "npc"), .BY[[1]], pos= 3, xpd= T)
      text(grconvertX(0.5, "npc"), grconvertY(1.32, "npc"), "x", pos= 3, xpd= T)
      text(grconvertX(0.75, "npc"), grconvertY(1.32, "npc"), .BY[[2]], pos= 3, xpd= T)
      pwm1 <- TF_clusters_PWMs$All_pwms_perc[[match(.BY[[1]], name(TF_clusters_PWMs$All_pwms_perc))]]
      my_seqlogo(as.matrix(pwm1), grconvertX(0.1, "npc"), grconvertY(1.05, "npc"), grconvertX(0.4, "npc"), grconvertY(1.35, "npc"))
      pwm2 <- TF_clusters_PWMs$All_pwms_perc[[match(.BY[[2]], name(TF_clusters_PWMs$All_pwms_perc))]]
      my_seqlogo(as.matrix(pwm2), grconvertX(0.6, "npc"), grconvertY(1.05, "npc"), grconvertX(0.9, "npc"), grconvertY(1.35, "npc"))
      
      .lm <- lm(diff~count1+count2, .SD)
      rsq <- paste0("R2", formatC(summary(.lm)$r.squared, format = "e", digits = 2))
      est1 <- paste0("mot1.est=", formatC(summary(.lm)$coefficients[1,1], format = "e", digits = 2))
      est2 <- paste0("mot2.est=", formatC(summary(.lm)$coefficients[2,1], format = "e", digits = 2))
      lines(box$stats[3,], lty= 2)
      legend("topleft", c(rsq, est1, est2), bty= "n")

      # Controls
      par(mar= c(3.5,5,0.5,2))
      my_boxplot(log2FoldChange~count1+count2, .SD, col_box = Cc)
      my_boxplot(median_L~count1+count2, .SD, col_box = Cc)
      my_boxplot(median_R~count1+count2, .SD, col_box = Cc)
      plot.new()
      print(.BY)
    }, .(mot1_s, mot2_s, effect)]
dev.off()