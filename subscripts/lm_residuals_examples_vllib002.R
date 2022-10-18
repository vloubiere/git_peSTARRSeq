setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- unique(dat[actClassL!= "inactive", .(L, indL)])[indL<=quantile(indL, 0.9), L]
QR <- unique(dat[actClassR!= "inactive", .(R, indR)])[indR<=quantile(indR, 0.9), R]
dat <- dat[L %in% QL & R %in% QR]
enh <- unique(rbind(dat[, .(ID= L)], dat[, .(ID= R)]))

# Select motifs and colors
# GATA/1; GATA/2; AP/1; Trl/4; DRE/1; DRE/2
sel <- c("cisbp__M5214", 
         "cisbp__M4690", 
         "stark__TGANTCA", 
         "stark__RSWGAGMRHRR", 
         "homer__AVYTATCGATAD_DREF", 
         "flyfactorsurvey__Dref_FlyReg_FBgn0015664")
col <- c("darkred", "red", "orange", "pink", "cyan", "blue")
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
counts <- vl_motif_counts(lib[enh, enh_sequence, on= "ID_vl==ID"], sel= sel)
enh <- cbind(enh, counts)
# Split enhancers based on motifs
dev <- enh[cisbp__M5214+cisbp__M4690+stark__TGANTCA+stark__RSWGAGMRHRR>1
           & homer__AVYTATCGATAD_DREF+flyfactorsurvey__Dref_FlyReg_FBgn0015664==0, ID]
hk <- enh[cisbp__M5214+cisbp__M4690+stark__TGANTCA+stark__RSWGAGMRHRR==0
          & homer__AVYTATCGATAD_DREF+flyfactorsurvey__Dref_FlyReg_FBgn0015664>0, ID]

#-----------------------------------------------#
# Select candidates
#-----------------------------------------------#
# Dev left enhancer
cand <- dat[predicted>4.5 & abs(indL-indR)<0.5 & L %in% dev]
cand[, gp:= fcase(log2FoldChange-predicted < -1 & R %in% hk, "Sub",
                  abs(log2FoldChange-predicted)<.5 & R %in% dev, "P",
                  log2FoldChange-predicted > 1 & R %in% dev, "Sup")]
cand <- cand[!is.na(gp)]
cand[, cut:= cut(predicted, seq(min(predicted), max(predicted), 0.25))] # Similar predicted
selL <- cand[, all(c("Sub", "P", "Sup") %in% gp), .(L, cut)][(V1), .(L, cut)]
cand <- cand[selL, on= c("L", "cut")]
cand[, gp:= factor(gp, c("Sub", "P", "Sup"))]
setorderv(cand, c("L", "cut", "gp"))

pdf("pdf/draft/lm_residuals_candidates_vllib002.pdf",
    height = 1.1*nrow(cand)/3)
par(mfrow= c(nrow(cand)/3, 3),
    oma= c(0,0,2.5,0),
    mar= c(3,10,1.5,0.5),
    tcl= -0.2,
    mgp= c(2,0.5,0),
    xpd= NA,
    font.main= 1,
    cex.main= 0.8)
cand[, {
  bar <- barplot(c(log2FoldChange, predicted, indR, indL),
                 horiz= T,
                 border= NA)
  mtext(paste0(gp, ": ", L, " x ", R),
        adj= 1,
        cex= 0.5)
  for(i in 1:2)
  {
    vl_plot_enh_motifs(vl_toDTranges(c(coorL, coorR)[i]),
                       genome= "dm3",
                       sel= sel,
                       xleft= par("usr")[1]-strwidth("M")*11,
                       ybottom= bar[c(4,3)[i],1]-0.5,
                       width= strwidth("M")*10,
                       height= 1,
                       p.cutoff = 5e-4,
                       col_alpha = 0.6,
                       col = col)
  }
  if(.GRP == 1)
    legend(grconvertX(0, "ndc", "user"),
           grconvertY(1, "ndc", "user")-strheight("M"),
           legend= vl_Dmel_motifs_DB_full[sel, motif_cluster, on= "motif_ID"],
           fill= col,
           xpd= NA,
           bty= "n",
           horiz= T,
           cex= 0.5)
  print("")
}, .(L, R, gp)]
dev.off()

#-----------------------------------------------#
# Plot example (see bottom page for the script used for selection)
#-----------------------------------------------#
# Select examples
ex <- dat[list("dev_medium_C_00549",
               c("control_flat_genomic_region_C_00780", "dev_weak_C_00395", "dev_weak_C_00268")), on= c("L", "R")]
ex[, Cc:=  c("cornflowerblue", "limegreen", "tomato")]
subSel <- sel
subCol <- col

pdf("pdf/draft/lm_residuals_examples_vllib002.pdf",
    height = 2.8,
    width = 4)
layout(matrix(c(1,1,2,3), nrow= 2), heights = c(0.3, 1))
par(mar= c(1,5,1,2),
    tcl= -0.2,
    las= 1,
    lwd= 0.5)
vl_boxplot(dat[L==ex$L[1], residuals],
           ylab= "Observed/Expected (log2)", 
           outline= T)
points(x= rep(1, 3), y= ex[, residuals])
par(mar= c(0,0,0,0),
    cex.axis= 0.7,
    tcl= -0.1,
    las= 1,
    xpd= NA,
    lwd= 0.5)
plot.new()
x <- 0.1
y <- seq(1-strwidth("M")*0.75,
         0+strwidth("M")*0.75, 
         length.out= length(subSel))
text(x, 
     y-0.01,
     vl_Dmel_motifs_DB_full[subSel, motif_cluster, on= "motif_ID"], 
     pos= 4,
     cex= 0.5)
points(rep(x, length(y)), 
       y, 
       pch= 15,
       xpd= NA,
       cex= 1,
       col= subCol)
vl_seqlogo(lapply(vl_Dmel_motifs_DB_full[subSel, pwms_perc, on= "motif_ID"], as.matrix),
           x= 0.275, 
           y = y, 
           pos = 4, 
           cex.width = 0.4,
           cex.height = 0.5)
par(mar= c(2,7,0.5,0.5))
ex[, {
  par(mgp= c(0.5, 0.5, 0))
  bar <- barplot(c(log2FoldChange[3],
                   predicted[3],
                   indR[3],
                   log2FoldChange[2],
                   predicted[2],
                   indR[2],
                   log2FoldChange[1],
                   predicted[1],
                   indR[1],
                   indL[1]),
                 col= c(rev(rep(Cc, each= 3)), "gold"),
                 names.arg= c(rep(c("5'/3' Obs. act.", "5'/3' Exp. lm", "3'"), 3), "5'"),
                 xlab= "Activity (log2)",
                 lwd= 0.1,
                 xaxt= "n",
                 xlim= c(0, 7),
                 horiz= T,
                 space= rep(c(1.5,0.25,0.25), length.out= 10),
                 xpd= NA)
  vl_plot_enh_motifs(bed = as.data.table(GRanges(c(coorL[1], coorR))),
                     genome= "dm3",
                     sel= subSel,
                     xleft= grconvertX(0, "nfc", "user"),
                     ybottom= bar[c(10,9,6,3),1]-0.5,
                     width= (par("usr")[1]-strwidth("   3'"))-grconvertX(0, "nfc", "user"),
                     height= 1,
                     p.cutoff = 5e-4,
                     col = subCol)
  par(mgp= c(0.75, 0, 0.15))
  axis(1, c(0,7), c(0,7))
  adj <- strwidth("M")*0.25
  for(i in 1:3)
    lines(c(log2FoldChange[i]+adj, 
            rep(max(c(log2FoldChange[i],predicted[i]))+adj*2, 2),
            predicted[i]+adj),
          rep(bar[c(c(7, 4, 1)[i], c(8, 5, 2)[i]), 1], each= 2))
}]
dev.off()
