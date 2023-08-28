dat <- readRDS("db/FC_tables/vllib006_DESeq2.rds")
dat[, multiplicative:= indL+indR]
dat[, additive:= log2(2^indL+2^indR)]
dat[, residuals:= log2FoldChange-multiplicative]

short <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
dat <- merge(dat,
             short[, .(L, R, indL, indR, log2FoldChange)],
             by= c("L", "R"),
             suffixes= c("_long", "_short"))

pdf("pdf/draft/Compare_add_mult_longSpacer_vllib006.pdf",
    width = 5,
    height = 2)
layout(matrix(c(1,2,3,3,4,4), nrow= 2),
       widths= c(0.45,1,1))
par(bty= "n",
    tcl= -0.1,
    las= 1,
    mar= c(3,3,0.2,0.2),
    oma= c(0,0,2,0),
    font.main= 1,
    mgp= c(0.75,0.25,0),
    cex= 0.5,
    cex.axis= 0.5,
    cex.lab= 0.6,
    cex.main= 0.8)
# Left
unique(dat[, .(L, indL_short, indL_long)])[,{
  plot(indL_short,
       indL_long,
       frame= F,
       pch= 16,
       col= adjustcolor("grey", 0.4),
       xaxt= "n",
       yaxt= "n",
       xlab= "Act. 300bp spacer (log2)",
       ylab= "Act. 2kb spacer (log2)")
  axis(1, lwd= 0.5)
  axis(2, lwd= 0.5)
  title(main= "5' cand.", xpd= NA, line= 0.25)
  PCC <- cor.test(indL_short,
                  indL_long)$estimate
  legend("topleft",
         legend= paste0("PCC= ", round(PCC, 2)),
         bty= "n",
         inset= c(-0.05, 0),
         cex= 0.6)
  abline(0, 1, lty= "11")
}]
# Right
unique(dat[, .(R, indR_short, indR_long)])[,{
  plot(indR_short,
       indR_long,
       frame= F,
       pch= 16,
       col= adjustcolor("grey", 0.4),
       xaxt= "n",
       yaxt= "n",
       xlab= "Act. 300bp spacer (log2)",
       ylab= "Act. 2kb spacer (log2)")
  axis(1, lwd= 0.5)
  axis(2, lwd= 0.5)
  title(main= "3' cand.", xpd= NA, line= 0.25)
  PCC <- cor.test(indR_short,
                  indR_long)$estimate
  legend("topleft",
         legend= paste0("PCC= ", round(PCC, 2)),
         bty= "n",
         inset= c(-0.05, 0),
         cex= 0.6)
  abline(0, 1, lty= "11")
}]
par(cex= 0.6,
    mgp= c(1.5,0.5,0),
    cex.axis= 1,
    cex.lab= 1,
    cex.main= 1)
# Pairs
dat[, {
  smoothScatter(additive,
                log2FoldChange_long,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Predicted additive (log2)",
                ylab= "Activity (log2)")
  title(main= "Additive", line= 1, xpd= NA)
  vl_plot_R2(rsquare = vl_model_eval(log2FoldChange_long, additive)$Rsquare,
             inset= c(-0.05, 0))
  clip(quantile(additive, 0.001),
       quantile(additive, 0.999),
       quantile(log2FoldChange_long, 0.001),
       quantile(log2FoldChange_long, 0.999))
  abline(0,1, lty= "11")
  smoothScatter(multiplicative,
                log2FoldChange_long,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Predicted multiplicative (log2)",
                ylab= "Activity (log2)")
  title(main= "Multiplicative", line= 1, xpd= NA)
  vl_plot_R2(rsquare = vl_model_eval(log2FoldChange_long, multiplicative)$Rsquare,
             inset= c(-0.05, 0))
  clip(quantile(multiplicative, 0.01),
       quantile(multiplicative, 0.99),
       quantile(log2FoldChange_long, 0.001),
       quantile(log2FoldChange_long, 0.999))
  abline(0,1, lty= "11")
}]
dev.off()