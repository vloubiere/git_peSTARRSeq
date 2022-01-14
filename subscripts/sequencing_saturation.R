setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

dat <- fread("Rdata/metadata_processed.txt")
dat <- dat[DESeq2 & file.exists(pairs_counts)]
dat <- dat[, fread(pairs_counts), .(group, pairs_counts, cdition, DESeq2_pseudo_rep)]

pdf("pdf/alignment/saturation.pdf", 4.5, 5)
par(las = 1)
dat[, {
  # Init
  plot(
    NA,
    xlim = c(0, 15),
    ylim = c(0, 0.6),
    xlab = "log2(counts+1)",
    ylab = "Density",
    main = group
  )
  # Legend
  leg <- .SD[, {
    paste0(cdition, "_rep", DESeq2_pseudo_rep)
  }, .(cdition, DESeq2_pseudo_rep)]$V1
  Cc <- ifelse(grepl("input", leg), "cornflowerblue", "tomato")
  Cc <- adjustcolor(Cc, 0.5)
  legend(
    "topright",
    legend = leg,
    col = Cc,
    lty = 1,
    lwd = 2,
    bty = "n"
  )
  # Denstiy replicates
  .SD[, {
    lines(density(log2(umi_counts + 1)),
          col = adjustcolor(Cc, 0.5)[.GRP],
          lwd = 2)
  }, .(cdition, DESeq2_pseudo_rep)]
  # Denstiy merge
  merge <- .SD[, .(umi_counts = sum(umi_counts)), .(cdition, L, R)]
  merge[, {
    lines(
      density(log2(umi_counts + 1)),
      col = ifelse(grepl("input", cdition), "cornflowerblue", "tomato"),
      lwd = 1,
      lty = 2
    )
  }, cdition]
  print(paste0(group, " -->>DONE!"))
}, group]
dev.off()