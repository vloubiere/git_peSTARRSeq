setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/vllib016_clustering_additive_scores_draft_figure.rds")

# Collapse per motif cluster
simp <- as.data.table(dat$mot$enr)
simp[vl_Dmel_motifs_DB_full, 
     c("motif_cluster_name", "motif_name", "pwm"):= .(i.Motif_cluster_name, i.motif_name, pwms_perc),
     on= "variable==uniqName_noSpecialChar"][]
sel <- simp[!is.na(log2OR), .(variable= variable[which.max(abs(log2OR))]), motif_cluster_name]
simp <- simp[variable %in% sel$variable & !is.na(log2OR)]

#----------#
act <- c(split(dat$rows$ind_act, dat$rows$cl), split(dat$cols$ind_act, dat$cols$cl))
names(act) <- paste0(rep(c("5' ", "3' "), each= 2), names(act)) 
#----------#
Cc <- circlize::colorRamp2(seq(-3,3, length.out= 20),
                           vl_palette_blueWhiteRed(20))

#----------------------------#
# PLOT
#----------------------------#
pdf("pdf/draft/Figure_3E.pdf",
    width= 5,
    height = 6)
par(las= 2,
    mar= c(3,17,0,5))
res <- plot.vl_enr_cl(simp,
                      padj_cutoff= 1e-5,
                      col= c("#00AEEF", "#EE2A7B"),
                      auto_margins= F, 
                      cex.balloons = 0.5)
logo <- as.data.table(res$x, keep.rownames= T)[, .(rn, bar= rev(.I))]
simp[logo, bar:= i.bar, on= "variable==rn"]
simp <- simp[!is.na(bar)]
simp[, top:= bar+0.35]
simp[, width:= strwidth("M")*1*ncol(as.matrix(pwm[[1]])), variable]
simp[, right:= par("usr")[1]-max(strwidth(paste0("MM", variable)))]
simp[, left:= right-width, variable]
simp[, height:= 0.7]
simp[, {
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
