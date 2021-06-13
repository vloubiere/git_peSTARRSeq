# Comparison Bernardo ####
lib <- readRDS("Rdata/uniq_library_final.rds")
res <- data.table(file= list.files("db/DE_analysis/", "add_scores.txt", full.names = T))
res[, lib:= tstrsplit(basename(file), "_", keep= 1)]
res <- res[, fread(file), res]
res[lib, single_L:= i.dev_log2FoldChange, on= "L==ID"]
res[lib, single_R:= i.dev_log2FoldChange, on= "R==ID"]

pdf("pdf/comparison_BA.pdf", width = 6, height = 10)
par(mfrow= c(3,2), pch= 16)
res[, {
  xl <- "twist-STARR-Seq"
  yl <- "pe-STARR-Seq"
  .c <- unique(.SD[, .(L, median_L, single_L)])
  .p <- seq(-2, 12, 0.1)
  plot(.c$median_L~.c$single_L, .c, main= paste0(lib, " L"), xlab= xl, ylab= yl)
  abline(0,1)
  lines(.p, predict(loess(.c$median_L~.c$single_L), .p), col= "red")
  .c <- unique(.SD[, .(R, median_R, single_R)])
  plot(median_R~single_R, .c, main= paste0(lib, " R"), xlab= xl, ylab= yl)
  lines(.p, predict(loess(.c$median_R~.c$single_R), .p), col= "red")
  abline(0,1)
}, lib]
dev.off()
####
