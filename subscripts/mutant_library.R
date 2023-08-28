setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/linear_models/FC_vllib029_mutLib_lm_predictions.rds")

# Extract data for combinations of interest ----
pl <- merge(dat[mutL=="WT" & mutR=="WT", .(IDL, IDR, indL, indR, log2FoldChange, residuals)],
            dat[mutL!="WT" & mutR!="WT", .(IDL, IDR, mutL, mutR, indL, indR, log2FoldChange, residuals)],
            by= c("IDL", "IDR"),
            suffixes= c("_wt", "_mut"))
pl[, cdition:= fcase(grepl("add.*Trl", mutL) & grepl("add.*Trl", mutR), "Added Trl motifs",
                     grepl("add.*Twist", mutL) & grepl("add.*Twist", mutR), "Added Twist motifs",
                     grepl("add.*Dref", mutL) & grepl("add.*Dref", mutR), "Added Dref motifs",
                     grepl("mut.*Trl", mutL) & grepl("mut.*Trl", mutR), "Mutated Trl motifs",
                     grepl("mut.*Twist", mutL) & grepl("mut.*Twist", mutR), "Mutated Twist motifs",
                     grepl("mut.*Dref", mutL) & grepl("mut.*Dref", mutR), "Mutated Dref motifs",
                     default= "mixed")]

pdf("pdf/draft/mutant_library_per_cdition.pdf", 5, 2.25)
layout(matrix(1:3, nrow = 1), 
       widths = c(0.75,1,0.5))
par(oma= c(0,0,4,6),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    las= 1,
    font.main= 1)
Cc <- c("lightgrey", "pink1")
pl[cdition!="mixed", {
  for(cutoff in c(-Inf, 0.5))
  {
    
    .SD[indL_wt>=cutoff & indR_wt>=cutoff, {
      par(mar= c(5,3.5,1.5,0.5))
      plot(residuals_wt,
           residuals_mut,
           xlab= "WT residuals (log2)",
           ylab= "Mut residuals (log2)",
           pch= 16,
           col= adjustcolor("grey", .2),
           bty= "n",
           cex= .5)
      PCC <- cor.test(residuals_wt,
                      residuals_mut)$estimate
      legend("topleft",
             legend= paste0("PCC= ", round(PCC, 2)),
             bty= "n",
             inset= c(-0.1, 0),
             cex= 0.6)
      clip(quantile(residuals_wt, 0.01),
           quantile(residuals_wt, 0.99),
           quantile(residuals_mut, 0.01),
           quantile(residuals_mut, 0.99))
      abline(0, 1, lty= "11")
      par(mar= c(2,3.5,0,0.5))
      vl_boxplot(indL_wt,
                 indL_mut,
                 indR_wt,
                 indR_mut,
                 log2FoldChange_wt,
                 log2FoldChange_mut,
                 compute_pval= list(c(1,2), c(3,4), c(5,6)),
                 ylab= "Activity (log2)",
                 notch= T,
                 xaxt= "n",
                 col= Cc)
      vl_boxplot(residuals_wt,
                 residuals_mut,
                 compute_pval= list(c(1,2)),
                 ylab= "Residuals (log2)",
                 xaxt= "n",
                 notch= T,
                 col= Cc)
      text(x= grconvertX(0.5, "ndc", "user"),
           y= grconvertY(1, "ndc", "user"),
           labels= paste0(unique(cdition),
                          ifelse(cutoff==(-Inf), " (all, n=", " (act, n="),
                          formatC(.N, big.mark = ","), ")"),
           pos= 1,
           xpd= NA,
           cex= 1.5)
      legend(par("usr")[2],
             par("usr")[4],
             fill= Cc,
             legend= c("WT", "Mut"),
             bty= "n",
             xpd= NA)
      .SD
    }]
  }
  .SD
}, cdition]
dev.off()

