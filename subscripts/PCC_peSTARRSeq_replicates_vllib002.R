setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]

pdf("pdf/draft/PCC_peSTARRSeq_replicates_vllib002.pdf", 
    height = 3, 
    width = 3)
par(las= 1,
    mar= c(3.5, 2.75, 0.5, 0.5),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    bty= "n")
dat[, {
  x <- log2(input_rep1+1)
  y <- log2(input_rep2+1)
  smoothScatter(log2(input_rep1+1), log2(input_rep2+1))
  legend("topleft",
         paste("PCC=", round(cor.test(x,y)$estimate, 2)),
         bty= "n")
  
  x <- log2(screen_rep1+1)
  y <- log2(screen_rep2+1)
  smoothScatter(log2(screen_rep1+1), log2(screen_rep2+1))
  legend("topleft", 
         paste("PCC=", round(cor.test(x,y)$estimate, 2)),
         bty= "n")
}]
dev.off()