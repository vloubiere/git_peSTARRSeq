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
dat <- dat[median_L>1 & median_R>1 & !is.na(diff) & !is.na(median_L) & !is.na(median_R)]
feat <- feat[ID %in% c(dat$enh_L, dat$enh_R)]
setkeyv(feat, "ID")
mot <- readRDS("Rdata/modeling_peSTARRSeq/motifs_linear_models_results_with_hk.rds")
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")

# Select motifs of interest
sel <- unique(melt(mot, measure.vars = patterns(mot= "^mot", Dmel= "^Dmel"))[, .(mot, Dmel)])
sel <- sel[data.table(Dmel= c("21.mor", "10.Dref_BEAF-32_grn_scro_pnr", "19.", "7.Ohler5(E-box)"),
                      side= rep(c("left", "right"), each=2), effect= rep(c("pos", "neg"), 2)), , on= "Dmel"]
# Add pe-STARR-Seq values and motif counts
sel <- sel[, dat[, .(enh_L, enh_R, diff, median_L, median_R)], (sel)]
setkeyv(sel, "side")
sel[side=="left", counts:= feat[enh_L][[col]], .(col= paste0("motif__", mot))]
sel[side=="right", counts:= feat[enh_R][[col]], .(col= paste0("motif__", mot))]
# "Cage"activity
sel <- sel[(side=="left" & between(median_L, 2, 3)) | (side=="right" & between(median_R, 2, 3))]

pdf("pdf/peSTARRSeq/examples_lm_motifs_LR_with_hk.pdf", 14, 5.5)
lay <- matrix(c(1,3,4,2,5,6,15,7,8,15,13,14,15,11,12,15,9,10), nrow= 3, byrow= T)
layout(lay, widths = c(0.4, 1, 1, 0.4, 1, 1), heights = c(1.3,1,1))
par(las= 1)

# All residuals 
par(mar= c(3.5,5,4,2))
sel[, 
    {
      .c <- .SD[, .(all= "all"), .(enh_L, enh_R, diff)]
      my_boxplot(diff~all , .c, outline= T, ylab= "Residuals (log2)", xlab= NA)
      mtext("all", 1, line= 1)
    }, side]

# Residuals
sel[, 
    {
      # Residuals
      my_boxplot(diff~counts, .SD, outline= T, ylab= "Residuals (log2)", xlab= NA, col_box = ifelse(effect=="pos", "tomato", "cornflowerblue"))
      .lm <- lm(diff~counts, .SD)
      abline(.lm, lty= 2)
      rsq <- paste("R2=", formatC(summary(.lm)$r.squared, format = "e", digits = 2))
      est <- paste("estimate=", formatC(summary(.lm)$coefficients[2,1], format = "e", digits = 2))
      legend("bottomright", legend = c(rsq, est), bty= "n")
      
      # Logo
      pwm <- TF_clusters_PWMs$All_pwms_perc[[match(mot, name(TF_clusters_PWMs$All_pwms_perc))]]
      my_seqlogo(as.matrix(pwm), grconvertX(0.4, "npc"), grconvertY(1.05, "npc"), grconvertX(1, "npc"), grconvertY(1.4, "npc"))
      text(grconvertX(0.4, "npc"), grconvertY(1.2, "npc"), Dmel, pos= 2, xpd= T)
    }, .(Dmel, mot, side, effect)]

# Individual activity constrained candidate
par(mar= c(3.5,5,0.5,2))
sel[,
    {
      if(side=="left")
      {
        .c <- unique(.SD[, .(enh_L, med= median_L, counts)])
        yl <- paste(side, "individual act. (log2)")
      }else if(side=="right")
      {
        .c <- unique(.SD[, .(enh_R, med= median_R, counts)])
        yl <- paste(side, "individual act. (log2)")
      }
      my_boxplot(med~counts, .c, ylim= c(0,8), ylab= yl, xlab= NA, col_box = ifelse(effect=="pos", "tomato", "cornflowerblue"))
      .lm <- lm(med~counts, .c)
      abline(.lm, lty= 2)
      rsq <- paste("R2=", formatC(summary(.lm)$r.squared, format = "e", digits = 2))
      est <- paste("estimate=", formatC(summary(.lm)$coefficients[2,1], format = "e", digits = 2))
      legend("topleft", legend = c(rsq, est), bty= "n")
      if(side=="right")
      {
        mtext("motif counts", side = 1, line= 2.2)
      }
      print("")
    }, .(side, Dmel, effect)]

# Individual activity paired candidate
par(mar= c(3.5,5,0.5,2))
sel[,
    {
      if(side=="left")
      {
        .c <- unique(.SD[, .(enh_R, med= median_R, counts)])
        yl <- "right individual act. (log2)"
      }else if(side=="right")
      {
        .c <- unique(.SD[, .(enh_L, med= median_L, counts)])
        yl <- "left individual act. (log2)"
      }
      my_boxplot(med~counts, .c, ylim= c(0,8), ylab= yl, xlab= NA, col_box = ifelse(effect=="pos", "tomato", "cornflowerblue"))
      .lm <- lm(med~counts, .c)
      abline(.lm, lty= 2)
      rsq <- paste("R2=", formatC(summary(.lm)$r.squared, format = "e", digits = 2))
      est <- paste("estimate=", formatC(summary(.lm)$coefficients[2,1], format = "e", digits = 2))
      legend("topleft", legend = c(rsq, est), bty= "n")
      if(side=="left")
      {
        mtext("motif counts", side = 1, line= 2.2)
      }
      print("")
    }, .(side, Dmel, effect)]
dev.off()
file.show("pdf/peSTARRSeq/examples_lm_motifs_LR_with_hk.pdf")
