setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)


pdf("pdf/modeling/barplot_motif_contribution.pdf", 
    height = 6)
par(mfrow= c(2, 1),
    mar= c(8,2,0.5,2),
    pty= "s",
    mgp = c(3, 0.5, 0))
mod <- readRDS("db/linear_models/global_models/vllib002_dev_x_dev__predict_log2FoldChange__actL*R+mergedMotSom.rds")
dat <- as.data.table(summary(mod)$coefficients, keep.rownames = T)
dat[, padj:= p.adjust(`Pr(>|t|)`, "fdr")]
dat[, Cc:= "lightgrey"]
dat[padj<0.05, Cc:= ifelse(`t value`>0, "tomato", "cornflowerblue")]
dat <- dat[order(padj)][1:10][order(`t value`)]
dat[, rn:= gsub("__", "/", gsub("motifSom__|_merge", "", rn))]
setorderv(dat, "t value")
barplot(dat$`t value`, 
        col= dat$Cc,
        names= dat$rn, 
        las= 2,
        cex.names= 0.5, 
        cex.axis= 0.5, 
        ylim= c(-100, 200), 
        axes= F)
axis(2, 
     las= 1, 
     line = 0.5, 
     tck= -0.025, 
     cex.axis= 0.6)
mtext("t-value", 
      2, 
      line = 2)

mod <- readRDS("db/linear_models/global_models/vllib016_hk_x_hk__predict_log2FoldChange__actL*R+mergedMotSom.rds")
dat <- as.data.table(summary(mod)$coefficients, keep.rownames = T)
dat[, padj:= p.adjust(`Pr(>|t|)`, "fdr")]
dat[, Cc:= "lightgrey"]
dat[padj<0.05, Cc:= ifelse(`t value`>0, "tomato", "cornflowerblue")]
dat <- dat[order(padj)][1:10][order(`t value`)]
dat[, rn:= gsub("__", "/", gsub("motifSom__|_merge", "", rn))]
dat[, rn:= gsub("Fer1_da/Fer2_da/Fer3_da/", "Fer/", rn)]
dat[, rn:= gsub("_da", "", rn)]
setorderv(dat, "t value")
barplot(dat$`t value`, 
        col= dat$Cc,
        names= dat$rn, 
        las= 2,
        cex.names= 0.5,
        cex.axis= 0.5, 
        ylim= c(-20, 30),
        axes= F)
axis(2, 
     las= 1, 
     line = 0.5, 
     tck= -0.025, 
     cex.axis= 0.6)
mtext("t-value", 
      2, 
      line = 2)
dev.off()

