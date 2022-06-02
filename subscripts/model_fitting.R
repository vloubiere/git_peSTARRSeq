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
clean <- vl_screen[vllib=="vllib002" & class== "enh./enh."]
clean <- clean[L %in% feat[group %in% c("hk", "dev", "shared"), ID] 
               & R %in% feat[group %in% c("hk", "dev", "shared"), ID]]



par(mfrow=c(4,2),
    las= 1,
    tcl= -0.2,
    mgp= c(2,0.5,0))
plot.new()
# Additive
smoothScatter(clean[, .("Expected additive= log2(2^5'+2^3')"= additive, 
                        "Observed"= log2FoldChange)],
              main= "Addtive model")
abline(0,1)

# Multiplicative
smoothScatter(clean[, .("Expected multiplicative= 5'+3'"= multiplicative, 
                        "Observed"= log2FoldChange)],
              main= "Multiplicative model")
abline(0,1)

# compare the two
smoothScatter(clean[, .("Observed-Expected additive"= log2FoldChange-additive,
                        "Observed-Expected multipalicative"= log2FoldChange-multiplicative)])
abline(h= 0, lty= "11")
abline(v= 0, lty= "11")
dens <- density(clean[, log2FoldChange-additive])
lines(dens$x, 
      par("usr")[4]+(dens$y/max(dens$y))*0.9*(grconvertY(1, "nfc", "user")-par("usr")[4]), xpd= T)
segments(median(clean[, log2FoldChange-additive]),
         par("usr")[4],
         median(clean[, log2FoldChange-additive]),
         par("usr")[4]+0.9*(grconvertY(1, "nfc", "user")-par("usr")[4]),
         lwd= 2,
         xpd= T)
dens <- density(clean[, log2FoldChange-multiplicative])
lines(par("usr")[2]+(dens$y/max(dens$y))*0.9*(grconvertX(1, "nfc", "user")-par("usr")[2]),
      dens$x, xpd= T)
segments(par("usr")[2],
         median(clean[, log2FoldChange-multiplicative]),
         par("usr")[2]+0.9*(grconvertX(1, "nfc", "user")-par("usr")[2]),
         median(clean[, log2FoldChange-multiplicative]),
         lwd= 2,
         xpd= T)





