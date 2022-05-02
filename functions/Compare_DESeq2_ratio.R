setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

dat <- fread("Rdata/metadata_processed.txt")
dat <- na.omit(unique(dat[, .(FC_file_DESeq, FC_file_ratio, group)]))
dat <- dat[file.exists(FC_file_DESeq)]
dat <- dat[, {
  merge(fread(FC_file_DESeq),
        fread(FC_file_ratio), 
        by= c("L", "R"), 
        suffix= c("_DESeq", "_ratio"))
}, (dat)]

pdf("pdf/alignment/DESeq2_vs_ratio_compare.pdf", width= 10, height = 3.5)
par(mfrow= c(1,3))
dat[, {
  smoothScatter(median_L_DESeq, 
                median_L_ratio,
                main= group)
  abline(0,1, lty= 2)
  smoothScatter(median_R_DESeq, 
                median_R_ratio,
                main= group)
  abline(0,1, lty= 2)
  smoothScatter(log2FoldChange_DESeq, 
                log2FoldChange_ratio,
                main= group)
  abline(0,1, lty= 2)
}, group]
dev.off()