setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

dev <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
dev <- dev[grepl("^control|^dev", L) & grepl("^control|^dev", R)]
prom <-  dev[grepl("^control", L) & grepl("^control|^dev", R)]
dev[, additive:= log2(2^indL+2^indR+2)]
dev[, additive:= log2(2^indL+2^indR-2^0)]
dev[, multiplicative:= log2(2^indL*2^indR/2^0)]

hk <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
hk <- hk[grepl("^control|^hk", L) & grepl("^control|^hk", R)]
hk[, additive:= log2(2^indL+2^indR-2^0)]
hk[, multiplicative:= log2(2^indL*2^indR/2^0)]

pdf("pdf/draft/scatterplot_dev_hk_models.pdf", 2.8, 3.2)
par(mfcol= c(2,2),
    cex= 9/12,
    # cex.axis= 6/8,
    mar= c(1.1, 1.1, 0, 0),
    oma= c(2, 2, 4, 0),
    mgp= c(1.25, 0.35, 0),
    tcl= -0.2,
    las= 1,
    bty= "n",
    cex.axis= 0.8)
for(mod in c("additive", "multiplicative"))
{
  for(class in c("dev", "hk"))
  {
    # Scatter plot
    x <- get(class)[[mod]]
    y <- get(class)[["log2FoldChange"]]
    col <- switch(class, "dev"= "springgreen", "hk"= "tomato")
    smoothScatter(x,
                  y,
                  colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                  col= "lightgrey",
                  xaxt= "n",
                  yaxt= "n")
    axis(1, padj= -0.6)
    axis(2)
    cx <- range(x)
    cy <- range(y)
    clip(cx[1], cx[2], cy[1], cy[2])
    abline(0, 1, lty= "11")
    
    # Add legend
    if(mod=="additive")
    {
      title(ylab= "Activity (log2)",
            xpd= NA)
    }
    if(class=="dev")
      title(main= mod,
            line= 1,
            xpd= NA)
    if(class=="hk")
      title(xlab= paste0("Predicted (log2)"),
            xpd= NA)
    
    # Rsquared
    rsqAddDev <- vl_model_eval(y, x)$Rsquare
    adj.rsqAdd <- 1-(((1-rsqAddDev)*(length(x)-1))/(length(x)-2-1))
    vl_plot_R2(rsquare = adj.rsqAdd,
               adjusted = F,
               cex= 0.8,
               inset= c(-0.125, 0))
  }
}
dev.off()