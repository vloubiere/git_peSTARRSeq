setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(vioplot)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")

# Import
dat <- lib[vllib=="vllib016"][, diff:= log2FoldChange-additive]
dat <- feat$add_feature(dat, feat$lib)
dat <- feat$add_feature(dat, feat$genomic)
dat <- feat$add_feature(dat, feat$top_motifs)
dat <- feat$add_feature(dat, feat$deepSTARR)
pred <- fread("/groups/stark/almeida/Projects/Vincent_pairedSTARRseq/results/gkmSVM_predictions/vl_pairs_sequences_predictions_hkCP.txt")
pred[, c("L", "R"):= tstrsplit(V1, "__")]
setnames(pred, "V2", "SVM_pred")
dat[pred, SVM:= SVM_pred, on= c("L", "R")]

# Filtering
L <- dat[group_L=="DHS_peak" & median_L<0
         & !group_R %in% c("control", "repressor") & act_wilcox_R<0.001]
R <- dat[!group_L %in% c("control", "repressor") & act_wilcox_L<0.001
         & group_R=="DHS_peak" & median_R<0]
LR <- dat[group_L=="DHS_peak" & median_L<0
          & group_R=="DHS_peak" & median_R<0]
LR[, dist:= log10(abs(end_L-start_R))]
LR[seqnames_L!=seqnames_R, dist:= max(LR$dist)]

# Pairs
matL <- dcast(L, median_L+L~R, value.var = "diff")[, L:= paste0(round(median_L, 1), "x ", L)][, !"median_L"]
matL <- as.matrix(matL, 1)
matL <- matL[, apply(matL, 2, function(x) sum(is.na(x))<5)]
matR <- dcast(R, median_R+R~L, value.var = "diff")[, R:= paste0(round(median_R, 1), "x ", R)][, !"median_R"]
matR <- as.matrix(transpose(matR, keep.names = "rn",  make.names = "R"), 1)
cmb <- as.matrix(dcast(LR, median_L+L~median_R+R, value.var = "diff")[, -1], 1)

pdf("pdf/analyses/RpS12_DHS_peaks.pdf", 
    width = 16, 
    height = 18)
layout(matrix(c(1,1,1,3,
                4,5,6,2,
                7,8,9,2,
                10,11,12,2), ncol= 4, byrow= T),
       widths= c(1,1,1,1.1),
       heights = c(0.8,1,1,1,1))
par(cex= 0.5,
    mar= c(12,12,5,0))
vl_heatmap(matL,
           cluster_rows= F,
           breaks = c(-2,0,2), 
           legend_title = "Obs./Add. (log2)",
           auto_margins = F)
par(mar= c(12,12,0,10))
vl_heatmap(matR,
           cluster_cols= F,
           breaks = c(-2,0,2), 
           legend_title = "Obs./Add. (log2)",
           auto_margins = F)
par(mar= c(12,12,5,10))
vl_heatmap(cmb,
           cluster_cols= F,
           cluster_rows= F,
           breaks = c(-2,0,2), 
           legend_title = "Obs./Add. (log2)",
           auto_margins = F,
           show_rownames = F,
           show_colnames = F)
par(mar= c(9,5,9,5),
    las= 1)
plot(L[, .(`median residuals inactDHS * active`= median(diff)), .(L, `deepSTARR prediction hk`= deep_hk_L)][, c(3,2)])
plot(R[, .(`median residuals active * inactDHS`= median(diff)), .(R, `deepSTARR prediction hk`= deep_hk_R)][, c(3,2)])
plot(L[, .(`median residuals inactDHS * active`= median(diff)), .(L, `deepSTARR prediction dev`= deep_dev_L)][, c(3,2)])
plot(R[, .(`median residuals active * inactDHS`= median(diff)), .(R, `deepSTARR prediction dev`= deep_dev_R)][, c(3,2)])
plot(LR[, .(`inactDHS x inactDHS pairs activity (log2)`= log2FoldChange, 
            `SVM prediction pair act.`= SVM)])
plot(LR[, .(`inactDHS x inactDHS pairs activity (log2)`= log2FoldChange, 
            `Genomic distance (log10)`= dist)])
dev.off()

test <- merge(dat[L=="DHS_peak_B_01067"], dat[R=="DHS_peak_B_01067"], by.x= "R", by.y= "L")
plot(test$log2FoldChange.x,
     test$log2FoldChange.y)
plot(test$diff.x,
     test$diff.y)
file.show("")