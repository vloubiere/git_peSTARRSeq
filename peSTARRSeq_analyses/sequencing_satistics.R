setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(gridExtra)

dat <- readRDS("db/Rdata/all_uniq_counts.rds")
dat <- dat[, .(UMI_collapsed_reads= sum(counts)), sample]
dat[, all_counts:= list.files("db/Rdata/", paste0(sample, ".*all.rds"), full.names = T), sample]
dat[, uniquely_aligned_reads:= nrow(readRDS(all_counts)), all_counts]

res <- melt(dat, id.vars= "sample", measure.vars = c("uniquely_aligned_reads", "UMI_collapsed_reads"))
res <- dcast(res, sample~variable)
res[, percentage_usable_reads:= round(UMI_collapsed_reads/uniquely_aligned_reads*100, 1)]
cols <- c("uniquely_aligned_reads", "UMI_collapsed_reads")
res[, (cols):= lapply(.SD, function(x) formatC(x, big.mark = ",")), .SDcols= cols]

pdf("pdf/peSTARRSeq/sequencing_statistics.pdf", 10, 5)
grid.table(res)
dev.off()