setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import Clustering object
cl <- readRDS("Rdata/vllib002_lm_residuals_hclust.rds")
rows <- cl$rows
cols <- cl$cols

# plot
pdf("pdf/draft/Figure_2B_hclust.pdf", 
    width= 7.6, 
    height = 7)
par(mar= c(3,3,3,7),
    mgp= c(3,0.1,0),
    tcl= -0.1,
    las= 1)

# Heatmap
plot(cl)

# Left individual activities
left <- grconvertX(0.6, "lines", "user")
width <- diff(c(left, par("usr")[1]-strwidth("M", cex= 0.6)))
rows[, rect_left:= (median-min(median))/diff(range(median))]
rows[, rect_left:= left+rect_left*width]
rows[, rect_right:= (0-min(median))/diff(range(median))]
rows[, rect_right:= left+rect_right*width]
rows[,{
  rect(rect_left,
       y-0.5,
       rect_right,
       y+0.5,
       border= NA, 
       xpd= T, 
       col= "grey40")
  ticks <- c(0, max(axisTicks(range(median), log= F)))
  at <- (ticks-min(median))/diff(range(median))
  at <- left+at*width
  axis(3, 
       at = at, 
       labels = ticks,
       xpd= T, 
       line= 0.25, 
       cex.axis= 0.5)
  segments(rect_right[1],
           par("usr")[3],
           rect_right[1],
           par("usr")[4],
           lwd= 0.1,
           xpd= T)
  text(mean(at),
       par("usr")[4],
       "median",
       xpd= T, 
       pos= 3, 
       offset= 1.25,
       cex= 0.7)
}]

# Right individual activities
par(mgp= c(3,0.25,0))
bottom <- grconvertY(0.6, "lines", "user")
height <- diff(c(bottom, par("usr")[3]-strheight("M", cex= 0.6)))
cols[, rect_top:= (median-min(median))/diff(range(median))]
cols[, rect_top:= bottom+rect_top*height]
cols[, rect_bot:= (0-min(median))/diff(range(median))]
cols[, rect_bot:= bottom+rect_bot*height]
cols[,{
  rect(x-0.5,
       rect_bot,
       x+0.5,
       rect_top,
       border= NA, 
       xpd= T, 
       col= "grey40")
  ticks <- c(0, max(axisTicks(range(median), log= F)))
  at <- (ticks-min(median))/diff(range(median))
  at <- bottom+at*height
  axis(4,
       at = at,
       labels = ticks,
       xpd= T,
       line= 0.25,
       cex.axis= 0.5)
  segments(par("usr")[1],
           rect_bot[1],
           par("usr")[2],
           rect_bot[1],
           lwd= 0.1,
           xpd= T)
  text(par("usr")[4],
       mean(at),
       "median",
       xpd= T,
       pos= 4,
       offset= 0.75,
       cex= 0.7)
}]
dev.off()