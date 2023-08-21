setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

dev <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
# dev <- dev[grepl("^control|^dev", L) & grepl("^control|^dev", R)]
dev <- dev[grepl("^dev", L) & grepl("^dev", R)]
dev[, additive:= log2(2^indL+2^indR)]
dev[, multiplicative:= indL+indR]

hk <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
# hk <- hk[grepl("^control|^hk", L) & grepl("^control|^hk", R)]
hk <- hk[grepl("^hk", L) & grepl("^hk", R)]
hk[, additive:= log2(2^indL+2^indR)]
hk[, multiplicative:= indL+indR]

pdf("pdf/draft/scatterplot_dev_hk_models.pdf", 2.8, 3.2)
par(mfcol= c(2,2),
    cex= 9/12,
    cex.axis= 6/8,
    mar= c(1.1, 1.1, 0, 0),
    oma= c(2, 2, 4, 0),
    mgp= c(1.25, 0.35, 0),
    tcl= -0.2,
    las= 1)
for(mod in c("additive", "multiplicative"))
{
  for(class in c("dev", "hk"))
  {
    # Scatter plot
    x <- get(class)[[mod]]
    y <- get(class)[["log2FoldChange"]]
    col <- switch(class, "dev"= "springgreen", "hk"= "tomato")
    vl_rasterScatterplot(x,
                         y,
                         col= adjustcolor(col, 0.2),
                         pch= 16)
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
par(mar= c(1.1, 3.1, 0, 1))
at <- c(1,2,4,5)
vl_boxplot(hk$log2FoldChange-hk$additive,
           hk$log2FoldChange-hk$multiplicative,
           dev$log2FoldChange-dev$additive,
           dev$log2FoldChange-dev$multiplicative,
           col= rep(c("tomato", "palegreen"), each= 2),
           xaxt= "n",
           at= at,
           ylab= "Residuals (log2)",
           compute_pval= list(c(1,2), c(3,4)),
           notch= T)
vl_tilt_xaxis(at,
              labels= rep(c("Additive", "Multiplicative"), 2))
legend(par("usr")[2],
       par("usr")[4],
       fill= c("tomato", "palegreen"),
       legend= c("Housekeeping", "Developmental"),
       xpd= NA,
       bty= "n")
abline(h= 0, lty= "11")
dev.off()