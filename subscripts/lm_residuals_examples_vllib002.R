setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & !grepl("^control", L) & !grepl("^control", R)]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
sequences <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
# Motifs to be plotted
# sel <- c("AP1/1__42", "Trl/1__1", "GATA/2__6", "Ebox/CATATG/twi__2", "Ebox/CACGTG/1__1", "DRE/1__2", "DRE/2")
# col <- c("darkred", "red", "tomato", "pink", "cyan", "cornflowerblue", "blue")
sel <- c("AP1/1__42", "Trl/1__1", "GATA/2__6", "DRE/1__2", "DRE/2")
col <- c("darkred", "red", "tomato", "cyan", "blue")
col <- adjustcolor(col, 0.5)
dat <- merge(lib,
             feat[, c("ID", sel), with= F],
             by.x= "L",
             by.y= "ID")
dat <- merge(dat,
             feat[, c("ID", sel), with= F],
             by.x= "R",
             by.y= "ID",
             suffixes= c("__L", "__R"))
dat[, predicted:= predict(model)]

#-----------------------------------------------#
# Select candidates
#-----------------------------------------------#
# cand <- dat[predicted>4.5
#             & abs(median_L-median_R)<0.5
#             & `AP1/1__42__L`+`Trl/1__1__L`+`GATA/2__6__L`+`Ebox/CATATG/twi__2__L`>0 # Dev motifs
#             & `DRE/1__2__L`+`DRE/2__L`+`Ebox/CACGTG/1__1__L`==0]
# cand[, gp:= fcase(log2FoldChange-predicted < -1 # Sub efficient
#                   & `DRE/1__2__R`+`DRE/2__R`+`Ebox/CACGTG/1__1__R`>1, "Sub",
#                   abs(log2FoldChange-predicted)<.5 # Well predicted
#                   & `DRE/1__2__R`+`DRE/2__R`+`Ebox/CACGTG/1__1__R`<2, "P",
#                   log2FoldChange-predicted > 1 # Synergystic
#                   & `AP1/1__42__R`+`Trl/1__1__R`+`GATA/2__6__R`+`Ebox/CATATG/twi__2__R`>0
#                   & `DRE/1__2__R`+`DRE/2__R`+`Ebox/CACGTG/1__1__R`==0, "Sup")]
# cand <- cand[!is.na(gp)]
# cand[, cut:= cut(predicted, seq(min(predicted), max(predicted), 0.25))]
# selL <- cand[, all(c("Sub", "P", "Sup") %in% gp), .(L, cut)][(V1), .(L, cut)]
# cand <- cand[selL, on= c("L", "cut")]
# cand[, gp:= factor(gp, c("Sub", "P", "Sup"))]
# setorderv(cand, c("L", "cut", "gp"))
# 
# pdf("pdf/draft/lm_residuals_examples_vllib002.pdf",
#     height = 1.1*nrow(cand)/3)
# par(mfrow= c(nrow(cand)/3, 3),
#     oma= c(0,0,2.5,0),
#     mar= c(3,10,1.5,0.5),
#     tcl= -0.2,
#     mgp= c(2,0.5,0),
#     xpd= NA,
#     font.main= 1,
#     cex.main= 0.8)
# cand[, {
#   seqL <- sequences[L, enh_sequence, on= "ID_vl"]
#   seqR <- sequences[R, enh_sequence, on= "ID_vl"]
#   bar <- barplot(c(log2FoldChange, predicted, median_R, median_L),
#                  horiz= T,
#                  border= NA)
#   mtext(paste0(gp, ": ", L, " x ", R),
#         adj= 1,
#         cex= 0.5)
#   for(i in 1:2)
#   {
#     vl_seqMotifs(c(seqL, seqR)[i],
#                  sel= sel,
#                  xleft= par("usr")[1]-strwidth("M")*11,
#                  ybottom= bar[c(4,3)[i],1]-0.5,
#                  width= strwidth("M")*10,
#                  height= 1,
#                  p.cutoff = 5e-4,
#                  col_alpha = 0.6,
#                  col = col)
#   }
#   if(.GRP == 1)
#     legend(grconvertX(0, "ndc", "user"),
#            grconvertY(1, "ndc", "user")-strheight("M"),
#            legend= sel,
#            fill= col,
#            xpd= NA,
#            bty= "n",
#            horiz= T,
#            cex= 0.5)
#   print("")
# }, .(L, R, gp)]
# dev.off()

#-----------------------------------------------#
# Plot example 
#-----------------------------------------------#
# Select examples
ex <- dat[list("dev_medium_B_00524",
               c("dev_weak_C_00365", "dev_strong_B_00302", "dev_medium_C_00542")), on= c("L", "R")]
ex[, Cc:=  c("cornflowerblue", "limegreen", "tomato")]

pdf("pdf/draft/lm_residuals_examples_vllib002.pdf",
    height = 3,
    width = 2.5)
layout(matrix(1:2, nrow= 2), heights = c(0.4, 1))
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
         length.out= length(sel))
text(x, 
     y-0.01,
     vl_Dmel_motifs_DB_full[sel, motif_cluster, on= "motif"], 
     pos= 4,
     offset= 0,
     cex= 0.7)
points(rep(x, length(y))-0.05, 
       y, 
       pch= 15,
       xpd= NA,
       cex= 1.5,
       col= col)
vl_seqlogo(lapply(vl_Dmel_motifs_DB_full[sel, pwms_perc, on= "motif"], as.matrix),
           x= 0.3, 
           y = y, 
           pos = 4, 
           cex.width = 0.8,
           cex.height = 0.9)
par(mar= c(2,7,0.5,0.5))
ex[, {
  par(mgp= c(0.5, 0.5, 0))
  bar <- barplot(c(log2FoldChange[3],
                   predicted[3],
                   median_R[3],
                   log2FoldChange[2],
                   predicted[2],
                   median_R[2],
                   log2FoldChange[1],
                   predicted[1],
                   median_R[1],
                   median_L[1]),
                 col= c(rep(Cc, each= 3), "gold"),
                 names.arg= c(rep(c("5'/3' Obs. act.", "5'/3' Exp. lm", "3'"), 3), "5'"),
                 xlab= "Activity",
                 border= NA,
                 xaxt= "n",
                 xlim= c(0, 8),
                 horiz= T,
                 space= rep(c(1.25,0.25,0.25), length.out= 10),
                 xpd= NA)
  vl_seqMotifs(sequences[c(L[1], R), enh_sequence, on= "ID_vl"],
               sel= sel,
               xleft= grconvertX(1, "line", "user"),
               ybottom= bar[c(10,9,6,3),1]-0.5,
               width= (par("usr")[1]-strwidth("   3'"))-grconvertX(1, "line", "user"),
               height= 1,
               p.cutoff = 5e-4,
               col = col)
  par(mgp= c(0.75, 0, 0.15))
  axis(1, c(0,8), c(0,8))
  adj <- strwidth("M")*0.25
  for(i in 1:3)
    lines(c(log2FoldChange[i]+adj, 
            rep(max(c(log2FoldChange[i],predicted[i]))+adj*2, 2),
            predicted[i]+adj),
          rep(bar[c(c(7, 4, 1)[i], c(8, 5, 2)[i]), 1], each= 2))
}]
dev.off()
