setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import dataset
if(!exists("feat"))
  feat <- readRDS("Rdata/final_300bp_enhancer_features_w_motifs.rds")
if(!exists("vl_screen"))
  vl_screen <- readRDS("Rdata/final_results_table.rds")

#-----------------------------------------------#
# Train activity based models and compute residuals
#-----------------------------------------------#
clean <- vl_screen[vllib=="vllib002"]

clean[, test:= rowMeans(data.frame((screen_rep1+1/sum(screen_rep1+1))/(input_rep1+1/sum(input_rep1+1)),
                                   (screen_rep2+1/sum(screen_rep2+1))/(input_rep2+1/sum(input_rep2+1))))]
clean[, test_L:= mean(test[grepl("^control", R)]), L]
clean[, test_R:= mean(test[grepl("^control", L)]), R]
smoothScatter(clean[class=="enh./enh.", .(log2(test_L+test_R), log2(test))])
plot(clean[class=="enh./enh."][sample(1e7, 10000), .(test_L+test_R, log(test))])
plot(clean[class=="enh./enh."][sample(1e7, 10000), .(test_L*test_R, log(test))])
abline(0, 1)

boxplot(clean[class=="enh./enh.", .(act= log2(test), add= log2(test_L+test_R), log2(test)-log2(test_L+test_R))], outline= F)
abline(h= 0, lty= "11")

par(mfrow=c(1,2))
smoothScatter(clean[class=="enh./enh.", .(log2FoldChange, log2(test))])
abline(0,1)
abline(-4.5,1)
smoothScatter(clean[class=="enh./enh.", .(log2(2^median_L+2^median_R), log2(test_L+test_R))])
abline(0,1)
abline(-4.5,1)


ex <- clean[L %in% c("dev_medium_A_00225", "control_Ecoli_B_00976") & R %in% c("dev_inactive_C_00134", "control_Ecoli_B_00991")]
barplot(as.matrix(ex[, .(L, input_rep1, screen_rep1)], 1), beside = T)
