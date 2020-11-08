setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
options(datatable.print.topn= 1)
require(data.table)
require(pheatmap)
require(circlize)

if(!file.exists("Rdata/modeling_peSTARRSeq/motifs_linear_models_results_with_hk.rds"))
{
  feat <- readRDS("Rdata/library/lib_features.rds")
  feat <- feat[(vl)]
  setkeyv(feat, "ID")
  
  dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
  dat <- dat[median_L>1 & median_R>1 & !is.na(diff) & !is.na(median_L) & !is.na(median_R)]
  
  #--------------------------------------------------------------#
  # lms Left~right motif combinations
  #--------------------------------------------------------------#
  cmbn <- CJ(mot1= grep("motif", colnames(feat), value= T),
             mot2= grep("motif", colnames(feat), value= T))
  cmbn[, c("mot1_s", "mot2_s") := lapply(.SD, function(x) gsub("^motif__", "", x))]
  som <- readRDS("Rdata/motifs/som_enriched_motifs.rds")$info
  cmbn[som, Dmel1:= i.Dmel_s, on= "mot1_s==best_match"]
  cmbn[som, Dmel2:= i.Dmel_s, on= "mot2_s==best_match"]
  cmbn <- cmbn <- cmbn[, .(.id= c("LR", "SUM"),
                           f1= "diff~median_L*median_R", 
                           f2= c("+counts1*counts2", "+countss1*countss2")), (cmbn)]
  cmbn[, median_L:= list(dat$median_L)]
  cmbn[, median_R:= list(dat$median_R)]
  cmbn[, diff:= list(dat$diff)]
  cmbn[, counts1:= .(list(feat[dat$enh_L][[mot1]])), mot1]
  cmbn[, counts2:= .(list(feat[dat$enh_R][[mot2]])), mot2]
  cmbn[, countss1:= .(list(feat[dat$enh_L][[mot1]]+feat[dat$enh_R][[mot1]])), mot1]
  cmbn[, countss2:= .(list(feat[dat$enh_L][[mot2]]+feat[dat$enh_R][[mot2]])), mot2]
  cmbn[, formula:= paste0(f1, f2), .(f1, f2)]
  
  res <- cmbn[, 
              {
                current <- as.data.table(lapply(.SD, unlist))
                data.table(summary(lm(as.formula(formula), current))$coefficients, keep.rownames = T)
              }, .(mot1_s, mot2_s, Dmel1, Dmel2, .id, formula), .SDcols= c("median_L", "median_R", "diff", "counts1", "counts2", "countss1", "countss2")]
  
  # SAVE
  saveRDS(res, "Rdata/modeling_peSTARRSeq/motifs_linear_models_results_with_hk.rds")
}
dat <- readRDS("Rdata/modeling_peSTARRSeq/motifs_linear_models_results_with_hk.rds")
dat[`Pr(>|t|)`>0.05, `t value`:=0]

#--------------------------------------------------------#
# PLOT
#--------------------------------------------------------#
Cc <- colorRamp2(c(-5, 0, 5), c("cornflowerblue", "white", "tomato"))

pdf("pdf/peSTARRSeq/motifs_lm_heatmap_with_hk.pdf", width = 15, height = 14)
layout(matrix(4:1, ncol=2))

par(xaxs= "i", yaxs= "i", mar= c(5,0.25,2,7))
.c <- dat[.id=="LR" & rn=="counts1:counts2"]
res <- my_heatmap(.c, "Dmel1", "Dmel2", value.var = "t value", row_labels = F, col_labels = F, leg_title = "LR motifs\ninteraction\n(t.stat)", 
                  leg_title_dist = 0.075, leg_y = 0.95, leg_x = 1.075, breaks = c(-5,0,5))
pval <- cut(res[["Pr(>|t|)"]], c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), c("****", "***", "**", "*", "N.S"))
text(res$xcoor, res$ycoor, pval, cex= ifelse(pval=="N.S", 0.5, 0.8))

# Right
par(mar= c(0.25,0.25,35,7))
.R <- dat[.id=="LR" & rn=="counts2"][res, coor:= i.xcoor, on= "Dmel2"]
my_boxplot(`t value`~coor, .R, las= 1, col_box = Cc(.R[, median(`t value`), keyby= coor]$V1), xaxt= "n", yaxt= "n", outline=T)
abline(h= 0, lty= 2)
axis(4, c(-10, 0, 10), labels= c(-10, 0, 15), las= 1)
mtext("Right candidate", 2, line = 0.25)
mtext("t value", 4, line = 2, las= 1)
text(seq(max(.R$coor)), grconvertY(1.05, from = "npc", "user"), .R[, unique(Dmel2), keyby= coor]$V1, srt= 60, pos= 4, xpd= T, offset= 0)

# Left
par(mar= c(5,37.5,2,0.25))
.L <- dat[.id=="LR" & rn=="counts1"][res, coor:= i.ycoor, on= "Dmel1"]
my_boxplot(`t value`~coor, .L, col_box = Cc(.L[, median(`t value`), keyby= coor]$V1), horizontal = T, xaxt= "n", yaxt= "n", xlab= "", outline=T)
abline(v= 0, lty= 2)
axis(1, c(-30, 0, 20), labels= c(-30, 0, 20), las= 1)
mtext("Left candidate", line = 0.25)
mtext("t value", 1, line = 2)
text(grconvertX(-0.05, from = "npc", "user"), seq(max(.L$coor)), .L[, unique(Dmel1), keyby= coor]$V1, pos= 2, xpd= T, offset= 0)

# Annot
plot.new()
text(0.95, -0.05, "t. stat", pos= 2, cex= 2, xpd= T, offset= 2)
dev.off()
