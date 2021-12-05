setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

dat <- fread("Rdata/metadata_processed.txt")
dat <- dat[vllib %in% c("vllib017", "vllib019", "vllib021") 
           & file.exists(pairs_counts) 
           & cdition=="input"]
dat[is.na(DESeq2_pseudo_rep), DESeq2_pseudo_rep:= 2]
dat <- dat[, fread(pairs_counts), .(vllib, pairs_counts, cdition, DESeq2_pseudo_rep)]

pdf("pdf/alignment/PCC_replicates_w_wo_intron.pdf")
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
                xlab= "intron removed",
                ylab= "intron not removed")
  legend("topleft", 
         paste0("PCC= ", round(cor.test(x, y)$estimate, 2)),
         bty= "n")
  print(title)
}, .(vllib, cdition)]
dev.off()
