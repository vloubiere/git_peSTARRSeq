setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]
dat[, basalMean:= mean(log2FoldChange[ctlL & ctlR])]
dat[, `Multiplicative model`:= log2(2^indL*2^indR/2^basalMean)]

pdf("pdf/draft/Compare_mult_linear_predictions_vllib002.pdf",
    height = 3, 
    width = 3)
par(mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
smoothScatter(dat[, .(`Multiplicative model`,
                      `Linear model`= predicted)],
              col= adjustcolor(blues9[9], .3),
              xaxt= "n")
axis(1, padj = -1.25)
vl_plot_coeff(value= cor.test(dat$`Multiplicative model`,
                              dat$predicted)$estimate,
              digits = 3,
              type = "pcc",
              inset= c(-.1, 0),
              cex= 7/12)
clip(-2,10,-2,10)
abline(0, 1, lty= "11")
dev.off()