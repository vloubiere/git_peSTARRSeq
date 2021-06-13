# multiple enhancer genes vs few enhancer genes ####
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
dat <- dat[!(spike_in)]
feat <- readRDS("Rdata/master_lib_features.rds")
dat[feat, gene_L:= me_symbol, on= "L==ID"]
dat[feat, gene_R:= me_symbol, on= "R==ID"]
pdf("pdf/multiple_enhancer_genes_vs_few_enhancer.pdf", height = 3.5)
par(mfrow= c(1, 3), las= 2, mar= c(7,4,2,2))
dat[, {
  close <- .SD[!is.na(gene_L) & gene_L==gene_R, log2FoldChange-add]
  far <- .SD[is.na(gene_L) | is.na(gene_R) | gene_L!=gene_R, log2FoldChange-add]
  boxplot(close, far, notch= T, names= c("mult. enh.", "few enh."), main= lib)
  text(0.5, max(close), labels = length(close), pos= 4, offset = 0)
  text(1.5, max(far), labels = length(far), pos= 4, offset = 0)
  print("")
}, lib]
dev.off()
####
