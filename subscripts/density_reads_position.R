setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(pracma)

# meta <- fread("Rdata/metadata_processed.txt")
# dat <- meta[, fread(umi_counts, nrow= 100000), .(group, cdition, umi_counts)]
# dat <- dat[!grepl("^ts", L) & !grepl("^ts", R)]

pdf("pdf/alignment/density_reads_positions.pdf")
par(mfrow=c(2,1))
dat[, {
  dens <- sapply(1:300, function(x) sum(pos_L==x))/.N
  barplot(dens, 
          ylim= c(-1,1),
          main= paste(cdition, group))
  text(x = mean(par("usr")[c(1, 2)]), 
       y= 1,
       labels = paste(which(dens>0.05), collapse = " "),
       pos= 1)
  dens <- sapply(1:300, function(x) sum(pos_R==x))/.N
  barplot(-dens, add=T)
  text(x = mean(par("usr")[c(1, 2)]), 
       y= -1,
       labels = paste(which(dens>0.05), collapse = " "),
       pos= 3)
  abline(h= 0.05)
  abline(h= -0.05)
}, keyby= .(group, cdition)]
dev.off()