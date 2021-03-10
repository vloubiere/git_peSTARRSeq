# Homotypic vs heterotypic pairs ####
dat <- data.table(file= list.files("db/DE_analysis/", "DE.txt", full.names = T))
dat <- dat[, fread(file), file]
pdf("pdf/homotypic_pairs.pdf", height = 3.5)
par(mfrow= c(1, 3), las= 2, mar= c(7,4,2,2))
dat[, {
  controls <- .SD[grepl("control", L) & grepl("control", R)]
  enh <- .SD[!grepl("control", L) & !grepl("control", R)]
  boxplot(controls[L==R, log2FoldChange], 
          controls[L!=R, log2FoldChange], 
          enh[L==R, log2FoldChange], 
          enh[L!=R, log2FoldChange], 
          main= basename(file), notch= T,
          names= c("homo_ctl", "hetero_ctl", "homo_enh", "hetero_enh"))
  pval <- wilcox.test(enh[L==R, log2FoldChange], enh[L!=R, log2FoldChange])$p.value
  text(3.5, 10, round(pval, 3))
  print("")
}, file]
dev.off()
####
