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
# Sample for CV
set.seed(1)
sel_L <- clean[, 1-.N/clean[,.N], L][, sample(L, round(.N/10), prob = V1)]
set.seed(1)
sel_R <- clean[, 1-.N/clean[,.N], R][, sample(R, round(.N/10), prob = V1)]
clean[, set:= ifelse(L %in% sel_L | R %in% sel_R, "test", "train")]
# Train linear model
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= clean[set=="train"])
rsq <- vl_model_eval(observed = clean[set=="test", log2FoldChange], 
                     predicted = predict(model, new= clean[set=="test"]))
saveRDS(model, "Rdata/CV_linear_model_vllib002.rds")
# Train non linear model
nls.model <- nls(log2FoldChange ~ SSasymp(multiplicative, Asym, R0, lrc), 
                 data = clean)
nls.rsq <- vl_model_eval(observed = clean[set=="test", log2FoldChange], 
                         predicted = predict(nls.model, new= clean[set=="test"]))
saveRDS(nls.model, "Rdata/CV_asymp_model_vllib002.rds")

# Extract linear model coefficients
eq <- summary(model)$coefficients[, "Estimate"]
eq <- round(eq, 2)
eq <- paste0("Activity= ", eq[1], " + ", eq[2], "*5' + ", eq[3], "*3' ", eq[4], "*5':3'")

pdf("pdf/draft/additve_vs_multiplicative.pdf", 
    height= 8,
    width= 8)
par(mfrow=c(2,2),
    las= 1,
    tcl= -0.2,
    mgp= c(2,0.5,0),
    mar= c(5,5,2,2))

# compare the two
smoothScatter(clean[, .("Observed-Expected additive"= log2FoldChange-additive,
                        "Observed-Expected multipalicative"= log2FoldChange-multiplicative)],
              xlim= c(-8, 6),
              ylim= c(-8, 6))
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

# Multiplicative
smoothScatter(clean[, .("Expected multiplicative= 5'+3'"= multiplicative, 
                        "Observed"= log2FoldChange)],
              main= "Multiplicative model")
abline(0,1)

# Additive
smoothScatter(clean[, .("Expected additive= log2(2^5'+2^3')"= additive, 
                        "Observed"= log2FoldChange)],
              main= "Additive model")
abline(0,1)

# Multiplicative model
smoothScatter(clean[, .("Predicted"= predict(model, newdata = clean), 
                        "Observed"= log2FoldChange)])
mtext(eq, line= 1.5)
abline(0,1)

# Non linear model fit
smoothScatter(clean[, .("Expected multiplicative= 5'+3'"= multiplicative, 
                        "Observed"= log2FoldChange)],
              main= "Multiplicative model")
abline(0,1)
lines(seq(0, 15, length.out= 100),
      predict(nls.model, newdata = data.frame(multiplicative= seq(0, 15, length.out= 100))),
      col= "red",
      lty= 2,
      lwd= 2)

# Non linear model prediction
coeffs <- round(summary(nls.model)$coeff[,"Estimate"], 1)
smoothScatter(clean[, .("Predicted"= predict(nls.model, newdata = clean), 
                        "Observed"= log2FoldChange)],
              main= paste0(coeffs["Asym"], "+(", coeffs["R0"], "-", coeffs["Asym"], ")*exp(-exp(", coeffs["lrc"], ")*multiplicative)"))
abline(0,1)
dev.off()


