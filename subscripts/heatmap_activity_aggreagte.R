setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
{
  dat <- readRDS("Rdata/final_results_table.rds")
  dat <- dat[L!=R]
}

# PLOT
pdf("pdf/aggregate_activity/aggregate_activity_additivity.pdf",
    height = 4.5,
    width = 14.5)
par(mfrow= c(1, 3))
dat[, {
  # Filter for minimum N
  gL <- .SD[, length(unique(L))>10, group_L][(V1), group_L]
  gR <- .SD[, length(unique(R))>10, group_R][(V1), group_R]
  .c <- .SD[group_L %in% gL & group_R %in% gR]
  # Activity
  .act <- dcast(data = .c, 
                formula= group_L~group_R, 
                value.var= "log2FC_merge",
                fun.aggregate= median, na.rm= T)
  .act <- as.matrix(.act, 1)
  # Diff
  .add <- dcast(data = .c, 
                formula= group_L~group_R, 
                value.var= "diff_merge",
                fun.aggregate= median, na.rm= T)
  .add <- as.matrix(.add, 1)
  # row_ord <- order(rowSums(.act))
  # col_ord <- order(colSums(.act))
  # .act <- .act[row_ord, ]
  # .act <- .act[, col_ord]
  # .add <- .add[row_ord, ]
  # .add <- .add[, col_ord]
  vl_heatmap(.act, 
             cluster_rows = F, 
             cluster_cols = F, 
             col = c("blue", "cornflowerblue", "yellow"),
             main= cdition, 
             display_numbers = T, 
             legend_title = "Activity (log2)")
  obj <- vl_heatmap(.add, 
                    cluster_rows = F, 
                    cluster_cols = F,
                    main= cdition, 
                    breaks = c(-1.5, -0.25, 0.25, 1.5),
                    col = c("cornflowerblue", "white", "white", "tomato"),
                    display_numbers = T, 
                    legend_title = "Additivity (log2)")
  # Balloon plots
  if(any(.act<0))
    .act <- .act-min(.act)
  .act <- t(.act)
  .add <- t(.add)
  plot.new()
  plot.window(xlim = c(1, nrow(.act)),
              ylim = c(ncol(.act), 1))
  title(cdition)
  segments(1, 
           seq(ncol(.act)), 
           nrow(.act),
           seq(ncol(.act)))
  segments(seq(nrow(.act)),
           1,
           seq(nrow(.act)),
           ncol(.act))
  Ccvec <- obj$col
  Ccvec[Ccvec=="white"] <- "lightgrey"
  Cc <- colorRamp2(obj$breaks, 
                   Ccvec)
  points(unlist(row(.act)),
         unlist(col(.act)),
         cex= unlist(.act), 
         xpd= T,
         pch= 19,
         col= Cc(unlist(.add)))
  axis(1, 
       at= seq(ncol(.act)),
       colnames(.act),
       las= 2,
       lwd= 0, 
       line= -0.25)
  axis(2, 
       at= seq(nrow(.act)),
       rownames(.act),
       las= 2,
       lwd= 0, 
       line= -0.25)
  # Legend color
  xleft <- grconvertX(1.05, "npc", "user")
  xright <- grconvertX(1.05+0.05, "npc", "user")
  ybottom <- grconvertY(0.7, "npc", "user")
  ytop <- grconvertY(0.95, "npc", "user")
  rasterImage(matrix(rev(Cc(seq(obj$breaks[1], obj$breaks[length(obj$breaks)], length.out = 101)))),
              xleft,
              ybottom,
              xright,
              ytop,
              xpd=T)
  ticks <- axisTicks(range(obj$breaks), log=F)
  ymin.ticks <- ybottom+(min(ticks)-obj$breaks[1])/diff(range(obj$breaks))*(ytop-ybottom)
  ymax.ticks <- ybottom+(max(ticks)-obj$breaks[1])/diff(range(obj$breaks))*(ytop-ybottom)
  text(xright,
       seq(ymin.ticks, ymax.ticks, length.out = length(ticks)),
       labels = ticks,
       pos=4,
       xpd= T,
       cex= 0.6,
       offset= 0.25)
  text(xleft,
       0.95+0.5*strheight("A"),
       labels = "Add. score (log2)",
       pos= 4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  # Legend balloons
  scale <- axisTicks(range(.act), log= F, nint = 4)
  maxBalloonInch <- strheight("A", units = "inches", cex= max(scale))*0.75
  bx <- mean(c(xleft, xright))
  btop <- grconvertY(0.5, "npc", "user")
  bbot <- grconvertY(0, "npc", "user")
  points(rep(bx, length(scale)),
         seq(bbot, btop, length.out = length(scale)), 
         xpd= T,
         col= "black",
         cex= scale,
         pch= 16)
  text(rep(bx, length(scale)),
       seq(bbot, btop, length.out = length(scale)),
       labels= scale,
       pos= 4,
       xpd= T,
       offset= 1.6,
       cex= 0.8)
  text(xleft,
       grconvertY(0.6, "npc", "user"),
       labels = "Activity (log2)",
       pos= 4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  print(cdition)
}, cdition]
dev.off()

