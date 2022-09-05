setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

dat <- fread("Rdata/metadata_processed.txt")
dat <- dat[DESeq2 & file.exists(pairs_counts)]
dat <- dat[, fread(pairs_counts), .(vllib, pairs_counts, cdition, DESeq2_pseudo_rep)]
dat <- dat[, check:= all(c(1, 2) %in% DESeq2_pseudo_rep), .(vllib, cdition)]
dat <- dat[(check)]

pdf("pdf/alignment/PCC_replicates.pdf", width = 10, height = 10)
par(mfrow= c(4,4),
    mgp= c(3,1,0),
    tcl= -0.2,
    mar= c(4,4,3,2))
dat[, {
  title <- paste(vllib, cdition)
  mat <- dcast(.SD, 
               L+R~DESeq2_pseudo_rep, 
               value.var= "umi_counts")
  x <- log2(mat[, `1`]+1)
  y <- log2(mat[, `2`]+1)
  smoothScatter(x, 
                y,
                main= title,
                las=1,
                xlab= "Replicate 1",
                ylab= "Replicate 2")
  legend("topleft", 
         paste0("PCC= ", round(cor.test(x, y)$estimate, 2)),
         bty= "n")
  print(title)
}, .(vllib, cdition)]
dev.off()
