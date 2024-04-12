setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
maxL <- max(dat$indL)
maxR <- max(dat$indR)
pl <- list("Weak 5' enh."= dat[between(indL, 1, 1.5), .(ref= L, refInd= indL, enh= R, ind= indR, log2FoldChange)],
           "Strong 5' enh."= dat[between(indL, maxL-1, maxL), .(ref= L, refInd= indL, enh= R, ind= indR, log2FoldChange)],
           "Weakest 3' enh."= dat[between(indR, 1, 1.5), .(ref= R, refInd= indR, enh= L, ind= indL, log2FoldChange)],
           "Strong 3' enh."= dat[between(indR, maxR-1, maxR), .(ref= R, refInd= indR, enh= L, ind= indL, log2FoldChange)])
pl <- rbindlist(pl, idcol = T)
pl[, n:= length(unique(ref)), .id]
pl[, min:= min(refInd), .id]
pl[, max:= max(refInd), .id]
pl[, name:= paste0(.id, "(", round(min, 1), "-", round(max, 1), ", n= ", n, ")"), .id]
pl <- pl[, .(log2FoldChange= mean(log2FoldChange)), .(name, n, min, max, enh, ind)]

pdf("pdf/draft/review_promoter_saturation.pdf",
    height = 3, 
    width = 12)
par(mfrow= c(1,4),
    mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
pl[, {
  plot(ind,
       log2FoldChange,
       xlab= paste0(ifelse(grepl("3'", name), "5'", "3'"), " activity (log2)"),
       xaxt= "n",
       ylab= "Mean combined activity (log2)",
       col= adjustcolor("tomato", 0.3),
       cex= 0.8,
       pch= 16,
       xlim= c(0.5, 7.5),
       ylim= c(0, 10.5),
       main= name)
  axis(1, padj= -1.25)
  rect(par("usr")[1],
       min,
       par("usr")[2],
       max,
       col= adjustcolor("lightgrey", .3),
       border= NA)
  ref.lab <- paste(ifelse(grepl("3'", name), "3'", "5'"), "act. range")
  text(par("usr")[2],
       mean(c(min, max)),
       ref.lab,
       pos= 2,
       cex= 7/12,
       offset= 0)
  # Loess
  .lo <- loess(log2FoldChange~ind, .SD[order(ind)])
  lines(.lo$x, .lo$fitted, lty= "11")
  print("")
}, .(name, min, max)]
dev.off()