setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
lib <- lib[(class_act_L=="inactive" & median_L<1) & (class_act_R=="inactive" & median_R<1)]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
lib[, pred:= predict(model, .SD)]
lib[, residuals:= log2FoldChange-pred]

matL <- dcast(lib,
              L~R, 
              value.var = "residuals")
matL <- as.matrix(matL, 1)
prcomp(matL)$x

lib[, marL:= sum(residuals), L]
lib[, marR:= sum(residuals), R]


mat <- dcast(lib,
             marL+L~marR+R, 
             value.var = "residuals")
mat <- as.matrix(mat[, -1], 1)
while(sum(is.na(mat))>0.05*nrow(mat)*ncol(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}

pdf("test/test2.pdf")
vl_heatmap(mat, 
           breaks= c(-2,0,2), 
           cluster_rows= F, 
           cluster_cols= F)
dev.off()

#-----------------------------------------------#
# Make object
#-----------------------------------------------#
dat <- merge(lib[class_act_L=="inactive" & median_L<1, .(res= mean(residuals)), .(name= L, median= median_L)],
             lib[class_act_R=="inactive" & median_R<1, .(res= mean(residuals)), .(name= R, median= median_R)],
             by= "name", 
             suffixes= c("L", "R"))
dat[, class:= fcase(resL>0 & resR>0, "Boost 5'&3'",
                    resL>0, "Boost 5'",
                    resR>0, "Boost 3'",
                    default="No boost")]
dat[, col:= switch(class, 
                   "Boost 5'&3'"= "limegreen",
                   "Boost 5'"= "red",
                   "Boost 3'"= "blue",
                   "No boost"= "grey"), class]
dat[, class:= factor(class, c("Boost 5'&3'", "Boost 5'", "Boost 3'", "No boost"))]
dat[, mean:= apply(.SD, 1, mean), .SDcols= c("medianL", "medianR")]
dat[, cex:= (mean-min(mean))/diff(range(mean))*1.5+0.3]

#-----------------------------------------------#
# Compute motif enrichment
#-----------------------------------------------#
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
sel <- names(feat)[names(feat) %in% vl_Dmel_motifs_DB_full$motif]
mot <- as.matrix(feat[dat, ..sel, on= "ID==name"])

#-----------------------------------------------#
# Plot
#-----------------------------------------------#
pdf("pdf/draft/median_residuals_inactive_elements.pdf", 3.5, 3)
par(mar= c(3,3,1,3),
    mgp= c(1.5,0.5,0),
    tcl= -0.2,
    las= 1)
dat[, plot(resR, 
           resL,
           pch= 19, 
           col= adjustcolor(col, 0.3),
           xlab= "inactive 3' mean(residuals) (log2)",
           ylab= "inactive 5' mean(residuals) (log2)",
           cex= cex)] 
vl_balloonskey(seq(min(dat$cex), max(dat$cex), length.out= 4),
               round(seq(min(dat$mean), max(dat$mean), length.out= 4), 1),
               left= par("usr")[2]+strwidth("M"),
               top= par("usr")[4]-strheight("M")*3.8, 
               main = "Mean\nind.\nact.\n", 
               main.cex = 0.7)
abline(h= 0, lty= "11")
abline(v= 0, lty= "11")
dev.off()

pdf("pdf/draft/motif_enrichment_inactive_elements.pdf", 6.5, 4)
par(mar= c(3,20,2,6),
    mgp= c(1.5,0.25,0),
    tcl= -0.2,
    las= 2)
enrL <- vl_motif_enrich(mot[dat$class %in% c("Boost 5'&3'", "Boost 5'"),],
                        mot[dat$class %in% c("No boost"),], 
                        collapse_clusters = vl_Dmel_motifs_DB_full[colnames(mot), motif_cluster, on="motif"],
                        padj_cutoff = 0.01, 
                        top_enrich = 10,
                        add_motifs = T, 
                        cex.height = 1.2, 
                        breaks = c(3, 8.5))
title(main= "Boost 5' vs. No boost")
enrR <- vl_motif_enrich(mot[dat$class %in% c("Boost 5'&3'", "Boost 3'"),],
                        mot[dat$class %in% c("No boost"),], 
                        collapse_clusters = vl_Dmel_motifs_DB_full[colnames(mot), motif_cluster, on="motif"],
                        padj_cutoff = 0.01,
                        top_enrich = 10,
                        add_motifs = T, 
                        cex.height = 1.2, 
                        breaks = c(3, 8.5))
title(main= "Boost 3' vs. No boost")
dev.off()

sequences <- readRDS("Rdata/vl_library_twist008_112019.rds")
sequences <- as.data.table(sequences)
sub <- lib[!is.na(resL) & !is.na(resR) & grepl("dev", L) & grepl("dev", R)]
leftSel <- sub[class_act_L=="inactive" & class_act_R=="inactive", any(log2FoldChange-pred>2), L][(V1), L]
sub <- sub[L %in% leftSel & class_act_R=="inactive"]
sub <- sub[order(2^log2FoldChange-2^pred, decreasing = T), .SD[c(1, .N)], L]
sel <- rbind(enrL, enrR)
sel <- sel[padj<1e-4]
sel <- sel[, .SD[which.max(abs(log2OR))], motif_ID]
Cc <- circlize::colorRamp2(c(-3,-0.1, 0.1, 3), c("royalblue", "cornflowerblue", "pink", "red"))
Cc <- adjustcolor(Cc(sel$log2OR), 0.5)

pdf("pdf/draft/examples_super_additive_inactive_pairs.pdf", 
    height = 1.1*nrow(sub)/2)
par(mfrow= c(nrow(sub)/2, 2),
    oma= c(0,0,2.5,0),
    mar= c(3,10,1.5,0.5),
    tcl= -0.2,
    mgp= c(2,0.5,0),
    xpd= NA,
    font.main= 1,
    cex.main= 0.8)
sub[, {
  seqL <- sequences[L, enh_sequence, on= "ID_vl"]
  seqR <- sequences[R, enh_sequence, on= "ID_vl"]
  bar <- barplot(c(log2FoldChange, pred, median_R, median_L), 
                 horiz= T, 
                 border= NA,
                 main= paste0(L, " x ", R))
  for(i in 1:2)
  {
    vl_seqMotifs(c(seqL, seqR)[i], 
                 sel= sel$motif_ID, 
                 xleft= par("usr")[1]-strwidth("M")*11,
                 ybottom= bar[c(4,3)[i],1]-0.5,
                 width= strwidth("M")*10,
                 height= 1, 
                 p.cutoff = 5e-4, 
                 col_alpha = 0.6, 
                 col = Cc)
  }
  if(.GRP == 1)
    legend(grconvertX(0, "ndc", "user"), 
           grconvertY(1, "ndc", "user")-strheight("M"),
           legend= sel$motif_ID,
           fill= Cc,
           xpd= NA,
           bty= "n", 
           horiz= T,
           cex= 0.5)
  print("")
}, .(L, R)]
dev.off()


