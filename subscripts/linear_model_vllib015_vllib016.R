setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

dev <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
dev <- dev[grepl("^control|^dev", L) & grepl("^control|^dev", R)]
dev[, additive:= log2(2^indL+2^indR)]
dev[, multiplicative:= indL+indR]

hk <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
hk <- hk[grepl("^control|^hk", L) & grepl("^control|^hk", R)]
hk[, additive:= log2(2^indL+2^indR)]
hk[, multiplicative:= indL+indR]

#-----------------------------------------------#
# Different models
#-----------------------------------------------#
rsqAddDev <- vl_model_eval(dev$log2FoldChange, dev$additive)$Rsquare
adj.rsqAdd <- 1-(((1-rsqAddDev)*(nrow(dev)-1))/(nrow(dev)-2-1))

rsqMultDev <- vl_model_eval(dev$log2FoldChange, dev$multiplicative)$Rsquare
adj.rsqMult <- 1-(((1-rsqMultDev)*(nrow(dev)-1))/(nrow(dev)-2-1))

rsqAddHk <- vl_model_eval(hk$log2FoldChange, hk$additive)$Rsquare
adj.rsqAdd <- 1-(((1-rsqAddHk)*(nrow(hk)-1))/(nrow(hk)-2-1))

rsqMultHk <- vl_model_eval(hk$log2FoldChange, hk$multiplicative)$Rsquare
adj.rsqMult <- 1-(((1-rsqMultHk)*(nrow(hk)-1))/(nrow(hk)-2-1))

pdf("pdf/draft/Compare_add_mult_vllib015_vllib016_smoothScatter.pdf",
    height = 3.75,
    width = 6)
par(font.main= 1,
    las= 1,
    lend= 2,
    cex= 8/12,
    cex.axis= 7/8,
    cex.main= 1,
    mfrow= c(2,3),
    mar= c(4.1,4.1,2,2),
    mgp= c(1.5, 0.325, 0),
    tcl= -0.2,
    bty= "n")
smoothScatter(dev$log2FoldChange,
              dev$additive,
              colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
              xlab= "Combined activity (log2)",
              ylab= "5' + 3' activities",
              main= "Additive")
abline(0, 1, lty= "11")
legend('bottomright', 
       legend= bquote(Adj.~R^2 == .(round(rsqAddDev, 2))),
       bty= "n",
       cex= 7/8)
smoothScatter(dev$log2FoldChange,
              dev$multiplicative,
              colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
              xlab= "Combined activity (log2)",
              ylab= "5' x 3' activities",
              main= "Multiplicative")
abline(0, 1, lty= "11")
legend('bottomright', 
       legend= bquote(Adj.~R^2 == .(round(rsqMultDev, 2))),
       bty= "n",
       cex= 7/8)
plot.new()
smoothScatter(hk$log2FoldChange,
              hk$additive,
              colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
              xlab= "Combined activity (log2)",
              ylab= "5' + 3' activities",
              main= "Additive")
abline(0, 1, lty= "11")
legend('bottomright', 
       legend= bquote(Adj.~R^2 == .(round(rsqAddHk, 2))),
       bty= "n",
       cex= 7/8)
smoothScatter(hk$log2FoldChange,
              hk$multiplicative,
              colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
              xlab= "Combined activity (log2)",
              ylab= "5' x 3' activities",
              main= "Multiplicative")
abline(0, 1, lty= "11")
legend('bottomright', 
       legend= bquote(Adj.~R^2 == .(round(rsqMultHk, 2))),
       bty= "n",
       cex= 7/8)
dev.off()