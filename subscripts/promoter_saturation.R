setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_actPairs_lm_predictions.rds")
pl <- list("Weakest 5' enhancer"= dat[L==dat[actL!="Inactive"][which.min(indL), L] & actR!="Inactive", .(ref= indL, ind= indR, log2FoldChange, predicted)],
           "Strongest 5' enhancer"= dat[L==dat[which.max(indL), L] & actR!="Inactive", .(ref= indL, ind= indR, log2FoldChange, predicted)],
           "Weakest 3' enhancer"= dat[actL!="Inactive" & R==dat[actR!="Inactive"][which.min(indR), R], .(ref= indR, ind= indL, log2FoldChange, predicted)],
           "Strongest 3' enhancer"= dat[actL!="Inactive" & R==dat[which.max(indR), R], .(ref= indR, ind= indL, log2FoldChange, predicted)])
pl <- rbindlist(pl, idcol = T)

pdf("pdf/draft/promoter_saturation.pdf", 3.5, 3.5)
par(mfrow=c(2,2),
    lend= 3,
    mar= c(3,4,3,0.5),
    cex= 8/12,
    cex.axis= 7/8,
    font.main= 1,
    cex.main= 1,
    las= 1,
    tcl= -0.2,
    bty= "n",
    mgp= c(1.5,0.35,0))
pl[, {
  plot(ind,
       log2FoldChange,
       col= adjustcolor("lightgrey", 0.5),
       cex= 0.8,
       pch= 16,
       xlim= c(0.5, 7.5),
       ylim= c(0, 10.5),
       main= .id,
       ylab= "Combined acitivity (log2)",
       xlab= paste0(ifelse(grepl("3'", .id), "5'", "3'"), " activity (log2)"))
  abline(h= ref, lwd= 0.5)
  ref.lab <- paste(ifelse(grepl("3'", .id), "3'", "5'"), "activity")
  text(par("usr")[2],
       ref-strheight("M"),
       ref.lab,
       pos= 2,
       cex= 0.9)
  # Loess
  .lo <- loess(log2FoldChange~ind, .SD[order(ind)])
  lines(.lo$x, .lo$fitted, lty= "11")
  print("")
}, .(.id, ref)]
dev.off()