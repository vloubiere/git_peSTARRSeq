setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data ----
dev <- readRDS("db/FC_tables/DSCP_focused_WT_DESeq2.rds")[grepl("^dev", L) & grepl("^dev", R)]
hk <- readRDS("db/FC_tables/RpS12_focused_WT_DESeq2.rds")[grepl("^hk", L) & grepl("^hk", R)]

# Individual activity value ----
value <- 3.5
set.seed(1)
selDev <- dev[round(indL, 1)==value & round(indR, 1)==value][sample(seq(.N), 1)]
set.seed(1)
selHk <- hk[round(indL, 1)==value & round(indR, 1)==value][sample(seq(.N), 1)]
selDev
selHk
sel <- rbindlist(list(dev= selDev,
                      hk= selHk),
                 idcol = "CP")
sel[, Additive:= log2(2^indL+2^indR-1)]
sel[, Multiplicative:= indL+indR]
setnames(sel,
         c("indL", "indR", "log2FoldChange"),
         c("5' enh.", "3' enh.", "Pair"))
sel <- melt(sel,
            id.vars = "CP",
            measure.vars = c("5' enh.", "3' enh.", "Pair", "Additive", "Multiplicative"))
setorderv(sel, c("CP", "variable"), c(-1, 1))
sel[, col:= c("tomato", "limegreen")[.GRP], CP]

addFC <- function(value,
                  bar,
                  obsIdx,
                  predIdx,
                  adj= c(1, 2, 1),
                  horiz= F)
{
  width <- if(horiz)
    strwidth("M")*0.25 else
      strheight("M")*0.25
  x <- rep(bar[c(obsIdx, predIdx)], each= 2)
  y <- c(value[obsIdx]+width*adj[1],
         max(value[c(predIdx, obsIdx)])+width*adj[2], 
         max(value[c(predIdx, obsIdx)])+width*adj[2], 
         value[predIdx]+width*adj[3])
  if(horiz)
  {
    y <- x
    x <- y
  }
  lines(x,
        y,
        xpd= T,
        lwd= .75)
  x <- mean(bar[c(obsIdx, predIdx)])
  y <- max(value[c(predIdx, obsIdx)])+width*adj[2]
  if(horiz)
  {
    y <- x
    x <- y
  }
  text(x,
       y,
       pos= ifelse(horiz, 4, 3),
       paste0("x", round(value[predIdx]/value[obsIdx], 2)),
       cex= 5/12,
       offset= .1,
       xpd= T)
}

# Plot ----
pdf("pdf/draft/review_examples_add_mult_hkCP_dCP.pdf", 2.9, 2.7)
vl_par(mgp= c(.75, .2, 0),
       cex.axis= 6/12)
sel[,{
  bar <- barplot(value,
                 col= adjustcolor(col, .5),
                 space= c(rep(0.1, 5), 1, rep(0.1, 4)),
                 border= NA,
                 ylab= "Activity (log2)")
  vl_tilt_xaxis(bar,
                labels = variable)
  # Add FC
  addFC(value, bar, 3, 4)
  addFC(value, bar, 3, 5, adj= c(5, 2, 1))
  addFC(value, bar, 8, 9)
  addFC(value, bar, 8, 10, adj= c(5, 3, 1))
  .SD[]
}]
dev.off()