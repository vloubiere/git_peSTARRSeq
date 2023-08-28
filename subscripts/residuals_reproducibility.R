setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Actual FC table
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
# dat <- dat[actL!="Inactive" & actR!="Inactive"]

# Repliactes counts
counts <- DESeq2::counts(readRDS("db/dds/vllib002.dds"), normalized= T)
counts <- as.data.table(counts, keep.rownames= T)
counts[, c("L", "R"):= tstrsplit(rn, "__")]

# FoldChange
counts[, log2FoldChange_rep1:= log2(screen_rep1/input_rep1)]
counts[, log2FoldChange_rep2:= log2(screen_rep2/input_rep2)]

# Compute individual act
counts[, indL1:= mean(log2FoldChange_rep1), L]
counts[, indR1:= mean(log2FoldChange_rep1), R]
counts[, indL2:= mean(log2FoldChange_rep2), L]
counts[, indR2:= mean(log2FoldChange_rep2), R]

# Remove pairs for which combined or ind act could not be computed accurately
counts <- counts[dat, on= c("L", "R")]

# linear model
lm1 <- lm(log2FoldChange_rep1~indL1*indR1, counts)
lm2 <- lm(log2FoldChange_rep2~indL2*indR2, counts)
counts[, log2FoldChange_merge:= dat$log2FoldChange]
counts[, predicted_rep1:= predict(lm1)]
counts[, predicted_rep2:= predict(lm2)]
counts[, residuals_rep1:= log2FoldChange_rep1-predicted_rep1]
counts[, residuals_rep2:= log2FoldChange_rep2-predicted_rep2]

cmb <- data.table(c("predicted_rep1", "residuals_rep1"),
                  c("predicted_rep2", "residuals_rep2"))

pdf("pdf/draft/residuals_reproductibility.pdf",
    height= 8)
par(mfrow=c(3,3),
    tcl= -0.2,
    bty= "n",
    mgp= c(1.5, 0.5, 0),
    las= 1)
cmb[, {
  x <- counts[[V1]]
  y <- counts[[V2]]
  smoothScatter(x,
                y,
                xlab= V1[1],
                ylab= V2[1],
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                col= "lightgrey")
  PCC <- cor.test(x, y)$estimate
  legend("topleft",
         paste0("PCC= ", round(PCC, 3)),
         bty= "n")
  .SD
}, .(V1, V2)]
dev.off()