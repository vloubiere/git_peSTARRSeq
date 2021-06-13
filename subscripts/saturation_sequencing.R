# Check saturation ####
sat <- data.table(file= list.files("db/count", "_umi_counts.txt$", full.names = T))
sat[, sample:= tstrsplit(basename(file), "_umi", keep= 1)]
sat[, lib:= tstrsplit(sample, "_", keep= 1)]
sat <- sat[, fread(file), sat]
sat <- sat[(exist)]

pdf("pdf/sequencing_saturation.pdf", 5, 5)
sat[, {
  .c <- .SD[, .(dens= list(density(log2(count)+1))), sample]
  cc1 <- colorRampPalette(c("tomato", "darkred"))(.c[grepl("input", sample), .N])
  cc2 <- colorRampPalette(c("cornflowerblue", "royalblue2"))(.c[grepl("DSCP", sample), .N])
  .c[grepl("input", sample), Cc:= cc1[.GRP], sample]
  .c[grepl("DSCP", sample), Cc:= cc2[.GRP], sample]
  .x <- .c[, range(dens[[1]]$x), sample][, range(V1)]
  .y <- .c[, range(dens[[1]]$y), sample][, range(V1)]
  plot(NA, xlim= .x, ylim= .y, las= 1, xlab= "log2(counts+1)", ylab= "density", main= lib)
  legend("topright", fill= .c$Cc, legend = .c$sample, bty= "n")
  .c[, lines(dens[[1]], col= Cc), .(sample, Cc)]
  print("")
}, lib]
dev.off()  

####