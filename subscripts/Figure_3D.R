setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/vllib016_clustering_additive_scores_draft_figure.rds")
#----------#
enr_L <- dat$mot_enr_L$enr
enr_L[vl_Dmel_motifs_DB_full, 
    c("motif_cluster_name", "motif_name", "pwm"):= .(i.Motif_cluster_name, i.motif_name, pwms_perc),
    on= "variable==uniqName_noSpecialChar"]
sel <- enr_L[!is.na(log2OR), .(variable= variable[which.max(abs(log2OR))]), motif_cluster_name]
enr_L <- enr_L[variable %in% sel$variable]
setnames(enr_L, 
         c("variable", "motif_cluster_name"), 
         c("uniq_name", "variable"))
enr_L[variable=="Ebox/CACGTG/SREBP", variable:= "Ebox/SREBP"]
#----------#
enr_R <- dat$mot_enr_R$enr
enr_R[vl_Dmel_motifs_DB_full, 
      c("motif_cluster_name", "motif_name", "pwm"):= .(i.Motif_cluster_name, i.motif_name, pwms_perc),
      on= "variable==uniqName_noSpecialChar"]
sel <- enr_R[!is.na(log2OR), .(variable= variable[which.max(abs(log2OR))]), motif_cluster_name]
enr_R <- enr_R[variable %in% sel$variable]
setnames(enr_R, 
         c("variable", "motif_cluster_name"), 
         c("uniq_name", "variable"))
enr_R[, variable:= gsub("CACGTG/|/CAGCTG/CACCTG|CAGCTG/", "", variable)]
#----------#
act_L <- split(dat$x, 
               dat$rows$cl)
#----------#
act_R <- split(t(dat$x), 
               dat$cols$cl)
#----------#
Cc <- circlize::colorRamp2(seq(-3,3, length.out= 20),
                           vl_palette_blueWhiteRed(20))

#----------------------------#
# PLOT
#----------------------------#
pdf("pdf/draft/Figure_3D.pdf",
    width= 8.3,
    height = 12.5)
layout(matrix(1:2, 
              nrow= 2), 
       heights = c(1,4))
# Left enhancer
par(las= 2,
    mar= c(1.5,28.7,2,5.7))
vl_boxplot(act_L,
           violin= T, 
           violcol= adjustcolor(Cc(sapply(act_L, median, na.rm= T)), 0.7),
           ylab= "Observed/Exp. Add.")
mtext("5' enhancer", 
      cex= 2,
      las= 1, 
      line= 0.5)
abline(h= 0, 
       lty= 2)
# Enrichment
par(las= 2,
    mar= c(10,30,0,7))
res <- plot(enr_L,
            padj_cutoff= 0.001,
            N_top= 10,
            col= c("#00AEEF", "#EE2A7B"),
            auto_margins= F, 
            cex.balloons = 0.8)
logo <- as.data.table(res$x, keep.rownames= T)[, .(rn, bar= rev(.I))]
enr_L[logo, bar:= i.bar, on= "variable==rn"]
enr_L <- enr_L[!is.na(bar)]
enr_L[, top:= bar+0.35]
enr_L[, width:= strwidth("M")*1.35*ncol(as.matrix(pwm[[1]])), variable]
enr_L[, right:= par("usr")[1]-max(strwidth(paste0("MM", variable)))]
enr_L[, left:= right-width, variable]
enr_L[, height:= 0.7]
enr_L[, {
  vl_seqlogo(as.matrix(pwm[[1]]),
             xleft = left,
             ytop = top,
             width = width,
             height= height)
  segments(left,
           top-height,
           left+width,
           top-height,
           lwd= 0.25,
           xpd= T,
           col= "grey")
}, .(variable, top, left, width, height)]

# Right enhancer
par(las= 2,
    mar= c(1.5,30.7,2,5.7))
vl_boxplot(act_R,
           violin= T, 
           violcol= adjustcolor(Cc(sapply(act_R, median, na.rm= T)), 0.7),
           ylab= "Observed/Exp. Add.")
mtext("3' enhancer", 
      cex= 2,
      las= 1, 
      line= 0.5)
abline(h= 0, 
       lty= 2)
# Enrichment
par(mar= c(23,32,0,7))
res <- plot(enr_R,
            padj_cutoff= 0.001,
            N_top= 10,
            col= c("#00AEEF", "#EE2A7B"),
            auto_margins= F, 
            cex.balloons = 0.65)
logo <- as.data.table(res$x, keep.rownames= T)[, .(rn, bar= rev(.I))]
enr_R[logo, bar:= i.bar, on= "variable==rn"]
enr_R <- enr_R[!is.na(bar)]
enr_R[, top:= bar+0.35]
enr_R[, width:= strwidth("M")*1.35*ncol(as.matrix(pwm[[1]])), variable]
enr_R[, right:= par("usr")[1]-max(strwidth(paste0("MM", variable)))]
enr_R[, left:= right-width, variable]
enr_R[, height:= 0.7]
enr_R[, {
  vl_seqlogo(as.matrix(pwm[[1]]),
             xleft = left,
             ytop = top,
             width = width,
             height= height)
  segments(left,
           top-height,
           left+width,
           top-height,
           lwd= 0.25,
           xpd= T,
           col= "grey")
}, .(variable, top, left, width, height)]
dev.off()
