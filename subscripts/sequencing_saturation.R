setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

dat <-
  read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
dat <- as.data.table(dat)
cols <- colnames(dat)
dat[, (cols) := lapply(.SD, function(x)
  ifelse(x == "NA", NA, x)), .SDcols = cols]

dat <- dat[, .(file = list.files(
  "db/merged_counts/",
  paste0(
    DESeq2_group,
    "_",
    cdition,
    "_rep",
    DESeq2_pseudo_rep,
    "_merged.txt"
  ),
  full.names = T
)),
.(DESeq2_group, cdition, DESeq2_pseudo_rep)]
dat <- dat[, fread(file), (dat)]
dat <- dat[type != "switched"]

dir.create("pdf/alignment",
           showWarnings = F)

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
    main = DESeq2_group
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
  print(paste0(DESeq2_group, " -->>DONE!"))
}, DESeq2_group]
dev.off()
