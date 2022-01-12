setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
{
  dat <- readRDS("Rdata/final_results_table.rds")
  dat <- dat[L!=R]
  # Filter for minimum N
  dat[, check:= .SD[, rep(.N>100, .N), .(group_L, group_R)]$V1, vllib]
  dat <- dat[(check), !"check"]
}

# PLOT
pl <- merge(dat[vllib=="vllib015", .(L, R, log2FoldChange, diff)], 
            dat[vllib=="vllib016"],
            by= c("L", "R"))

pdf("pdf/wip.pdf", 
    width = 14, 
    height = 17)
par(mfrow= c(7, 7))
pl[, {
  plot(log2FoldChange.x,
       log2FoldChange.y,
       xlim= c(-5, 10),
       ylim= c(-5, 10),
       pch= 16,
       cex= 0.5,
       col= adjustcolor(colorRampPalette(c(col_L, col_R))(3)[2], 0.3),
       main= paste(group_L, group_R),
       xlab= "DSCP",
       ylab= "RpS12")
  abline(0, 1, lty= 2)
}, .(group_L, group_R, col_L, col_R)]
dev.off()

pdf("pdf/wip.pdf", 
    width = 14, 
    height = 17)
par(mfrow= c(7, 7))
pl[, {
  plot(log2((2^log2FoldChange.x+2^log2FoldChange.y)/2),
       log2FoldChange.x-log2FoldChange.y,
       xlim= c(-2, 8),
       ylim= c(-5, 5),
       pch= 16,
       cex= 0.5,
       col= adjustcolor(colorRampPalette(c(col_L, col_R))(3)[2], 0.3),
       main= paste(group_L, group_R),
       xlab= "o/e",
       ylab= "DSCP/RpS12")
  abline(h= 0, lty= 2)
}, .(group_L, group_R, col_L, col_R)]
dev.off()
