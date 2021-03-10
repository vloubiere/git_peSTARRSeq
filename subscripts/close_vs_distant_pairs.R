# Close versus distant ####
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
dat <- dat[!(spike_in)]
feat <- readRDS("Rdata/master_lib_features.rds")
dat[feat, c("seqnames_L", "pos_L", "gene_L"):= .(i.seqnames, round((i.end+i.start)/2), symbol), on= "L==ID"]
dat[feat, c("seqnames_R", "pos_R", "gene_R"):= .(i.seqnames, round((i.end+i.start)/2), symbol), on= "R==ID"]
dat[seqnames_L==seqnames_R, dist:= abs(pos_L-pos_R)]
dat[seqnames_L!=seqnames_R, dist:= Inf]
pdf("pdf/close_vs_distant_pairs.pdf", height = 3.5)
par(mfrow= c(1, 3), las= 2, mar= c(7,4,2,2))
cutoff <- 30000
dat[, {
  close <- .SD[L!=R & dist<cutoff, .(L, R, median_L, median_R, diff= log2FoldChange-add)]
  rest <- .SD[dist>=cutoff]
  far <- data.table()
  for(i in seq(nrow(close)))
  {
    sel <- which.min(abs(close[i, median_L]-rest$median_L)+abs(close[i, median_R]-rest$median_R))
    far <- rbind(far, rest[sel])
    rest <- rest[-sel]
  }
  close <- close$diff
  far <- far[, log2FoldChange-add]
  boxplot(close, far, notch= T, names= c("close", "far"), main= lib)
  text(0.5, max(close), labels = length(close), pos= 4, offset = 0)
  text(1.5, max(far), labels = length(far), pos= 4, offset = 0)
  print("")
}, lib]
dev.off()
####
