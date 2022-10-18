setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
pl <- rbindlist(list("Weakest active 5' enhancer"= dat[!ctlL & !ctlR & actClassL=="active"][L==L[which.min(indL)], .(ref= indL, ind= indR, log2FoldChange, predicted)],
                     "Strongest active 5' enhancer"= dat[!ctlL & !ctlR & actClassL=="active"][L==L[which.max(indL)], .(ref= indL, ind= indR, log2FoldChange, predicted)],
                     "Weakest active 3' enhancer"= dat[!ctlR & actClassR=="active"][R==R[which.min(indR)], .(ref= indR, ind= indL, log2FoldChange, predicted)],
                     "Strongest active 3' enhancer"= dat[!ctlR & actClassR=="active"][R==R[which.max(indR)], .(ref= indR, ind= indL, log2FoldChange, predicted)]), 
                idcol = T)

pdf("pdf/draft/promoter_saturation.pdf", 5, 5.5)
par(mfrow=c(2,2),
    mar= c(3,4,3,0.5),
    las= 1,
    tcl= -0.2,
    bty= "n",
    mgp= c(1.5,0.5,0))
pl[, {
  xl <- range(ind)
  yl <- range(c(log2FoldChange, predicted))
  smoothScatter(ind,
                log2FoldChange,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlim= xl+diff(xl)*c(-0.01, 0.01),
                ylim= yl+diff(yl)*c(-0.01, 0.01),
                main= .id,
                ylab= "Activity of the pair (log2)",
                xlab= paste0(ifelse(grepl("3'", .id), "5'", "3'"), " individual activity (log2)"))
  abline(h= ref[1], lty= 2)
  ref.lab <- paste0(ifelse(grepl("3'", .id), "3'", "5'"), " individual\nactivity")
  text(par("usr")[2]-strwidth(ref.lab)/2,
       ref[1],
       ref.lab,
       pos= 1,
       cex= 0.9)
  lines(ind, predicted, col= "red")
  ""
}, .id]
dev.off()