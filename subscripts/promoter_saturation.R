setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")
pl <- list("Weakest 5' enhancer"= dat[L==dat[actL!="Inactive"][which.min(indL), L] & actR!="Inactive", .(ref= indL, ind= indR, log2FoldChange, predicted)],
           "Strongest 5' enhancer"= dat[L==dat[which.max(indL), L] & actR!="Inactive", .(ref= indL, ind= indR, log2FoldChange, predicted)],
           "Weakest 3' enhancer"= dat[actL!="Inactive" & R==dat[actR!="Inactive"][which.min(indR), R], .(ref= indR, ind= indL, log2FoldChange, predicted)],
           "Strongest 3' enhancer"= dat[actL!="Inactive" & R==dat[which.max(indR), R], .(ref= indR, ind= indL, log2FoldChange, predicted)])
pl <- rbindlist(pl, idcol = T)

pdf("pdf/draft/promoter_saturation.pdf", 5, 5.5)
par(mfrow=c(2,2),
    mar= c(3,4,3,0.5),
    las= 1,
    tcl= -0.2,
    bty= "n",
    mgp= c(1.5,0.5,0))
pl[, {
  smoothScatter(ind,
                log2FoldChange,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlim= c(0.5, 7.5),
                ylim= c(0, 10.5),
                main= .id,
                ylab= "Activity of the pair (log2)",
                xlab= paste0(ifelse(grepl("3'", .id), "5'", "3'"), " individual activity (log2)"))
  abline(h= ref, lty= 2)
  ref.lab <- paste(ifelse(grepl("3'", .id), "3'", "5'"), "activity")
  text(par("usr")[2],
       ref-strheight("M"),
       ref.lab,
       pos= 2,
       cex= 0.9)
  # Loess
  .lo <- loess(log2FoldChange~ind, .SD[order(ind)])
  lines(.lo$x, .lo$fitted, col= "red")

  # # Linear model
  # .lm <- lm(log2FoldChange~ind)
  # abline(.lm, col= "red")
  # leg <- vl_model_equation(.lm, digits= 2)
  # leg <- gsub("log2FoldChange", "y", leg)
  # leg <- gsub("ind", "x", leg)
  # legend("topleft", leg, bty= "n")
  print("")
}, .(.id, ref)]
dev.off()

file.show("pdf/draft/promoter_saturation.pdf")