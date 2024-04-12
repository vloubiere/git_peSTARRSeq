setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
pl <- list("Weakest 5' enhancer"= dat[L==dat[actL!="Inactive"][which.min(indL), L] & actR!="Inactive", .(ref= indL, ind= indR, log2FoldChange, predicted)],
           "Strongest 5' enhancer"= dat[L==dat[which.max(indL), L] & actR!="Inactive", .(ref= indL, ind= indR, log2FoldChange, predicted)],
           "Weakest 3' enhancer"= dat[actL!="Inactive" & R==dat[actR!="Inactive"][which.min(indR), R], .(ref= indR, ind= indL, log2FoldChange, predicted)],
           "Strongest 3' enhancer"= dat[actL!="Inactive" & R==dat[which.max(indR), R], .(ref= indR, ind= indL, log2FoldChange, predicted)])
pl <- rbindlist(pl, idcol = T)

pdf("pdf/draft/promoter_saturation.pdf",
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
       xlab= paste0(ifelse(grepl("3'", .id), "5'", "3'"), " activity (log2)"),
       xaxt= "n",
       ylab= "Combined acitivity (log2)",
       col= adjustcolor("lightgrey", 0.5),
       cex= 0.8,
       pch= 16,
       xlim= c(0.5, 7.5),
       ylim= c(0, 10.5),
       main= .id)
  axis(1, padj= -1.25)
  abline(h= ref, lwd= 0.5)
  ref.lab <- paste(ifelse(grepl("3'", .id), "3'", "5'"), "activity")
  text(par("usr")[2],
       ref-strheight("M")/1.5,
       ref.lab,
       pos= 2,
       cex= 7/12)
  # Loess
  .lo <- loess(log2FoldChange~ind, .SD[order(ind)])
  lines(.lo$x, .lo$fitted, lty= "11")
  print("")
}, .(.id, ref)]
dev.off()