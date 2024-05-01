setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data with activity and residuals ----
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
dat$ctlL <- dat$ctlR <- NULL

# Add motif counts ----
mot <- readRDS("db/motif_counts/twist008_motif_counts_selected.rds")
dat <- merge(dat, mot, by.x= "L", by.y= "ID", all.x= T, allow.cartesian = T)
dat <- merge(dat, mot, by.x= "R", by.y= "ID", all.x= T, allow.cartesian = T, suffixes= c("__L", "__R"))

# Compute residuals/activity for each motif combination and t.test pval ----
cmb <- CJ(grep("__L$", names(dat), value = T),
          grep("__R$", names(dat), value = T))
res <- cmb[, {
  if(.GRP %% 100==0)
    print(.GRP)
  # Check motifs
  mot <- dat[[V1]]>0 & dat[[V2]]>0 # Expected motifs on both sides
  act_mot <- dat[(mot), log2FoldChange]
  res_mot <- dat[(mot), residuals]
  noMot <- dat[[V1]]==0 & dat[[V2]]==0 # No motifs
  act_nomot <- dat[(noMot), log2FoldChange]
  res_nomot <- dat[(noMot), residuals]
  # Return values
  .(act_mot= mean(act_mot),
    act_nomot= mean(act_nomot),
    act_pval= t.test(act_mot, act_nomot)$p.value,
    res_mot= mean(res_mot),
    res_nomot= mean(res_nomot),
    res_pval= t.test(res_mot, res_nomot)$p.value)
}, .(V1, V2)]
# FDR ----
res[, act_padj:= p.adjust(act_pval, "fdr")]
res[, res_padj:= p.adjust(res_pval, "fdr")]
plot(density(-log10(res$res_padj)),
     xlab= "-log10(padj)",
     col= "red",
     frame= F)
lines(density(-log10(res$act_padj)))
legend("topright",
       col= c("red", "black"),
       legend= c("Residuals", "Activity"),
       lwd= 1,
       bty= "n")
# Simplify names ----
res[, V1:= gsub("/", ".", V1)]
res[, V2:= gsub("/", ".", V2)]
# Compute difference between pairs with motifs and without motifs ----
res[, `Mean activity [motif / no motif] (log2)`:= act_mot-act_nomot]
res[, `Mean residuals [motif / no motif] (log2)`:= res_mot-res_nomot]

saveRDS(res, "db/motif_effect/motif_impact_act_res.rds")
