setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(diagram)
require(bezier)
require(plotrix)
require(vlfunctions)

pdf("pdf/draft/Sketch_peSTARRSeq.pdf", width = 4, height = 4)
par(mar= c(0,0,0,0), 
    xaxt= "n", 
    yaxt= "n", 
    lend= 3, 
    lwd= 1.3,
    cex= 1)

# Plot init
plot.new()

# Backbone
segments(0, 0.8, 0.3, 0.8)
segments(0.7, 0.8, 1, 0.8)
segments(0.3, 0.8, 0.335, 0.8, col= "red")
segments(0.665, 0.8, 0.7, 0.8, col= "blue")

# Promoter
segments(0.1, 0.8, 0.1, 0.85)
segments(0.1, 0.85, 0.15, 0.85, lend=2)
polygon(c(0.15,0.175,0.15),
        c(0.86,0.85,0.84), 
        border= NA, 
        col= "black")

# Boxes
rect(0.05, 0.775, 0.15, 0.825, col= "white")
text(0.1, 0.8, "CP")
rect(0.75, 0.775, 0.9, 0.825, col= "white")
text(0.825, 0.8, "pA site")

# Enhancer pairs
Cc <- adjustcolor(vl_palette_few_categ(7), 0.7)
pos <- seq(0, 0.1, length.out= length(Cc))
for(i in seq(pos))
{
  segments(0.35, c(0.8-pos[i]), 0.45, c(0.8-pos[i]), col= Cc[i], lwd= 3)
  segments(0.45, c(0.8-pos[i]), 0.55, c(0.8-pos[i]))
  segments(0.55, c(0.8-pos[i]), 0.65, c(0.8-pos[i]), col= rev(Cc)[i], lwd= 3)
}
text(0.4, 0.805, "5'seq", pos= 3, offset= 0.25)
text(0.6, 0.805, "3'seq", pos= 3, offset= 0.25)
text(0.5, 0.805, "spacer", pos= 3, offset= 0.25, cex= 0.7)

# Plasmids
px <- cos(seq(0, 2*pi, length.out= 100))
px <- px-min(px)
px <- px/max(px)/20
py <- sin(seq(0, 2*pi, length.out= 100))
py <- py-min(py)
py <- py/max(py)/20
set.seed(1)
ppx <- jitter(c(0.1,0.175, 0.25, 0.1375, 0.2125, 0.1,0.175))-0.1
set.seed(8)
ppy <- jitter(c(0.6,0.6,0.6,0.65,0.65,0.7,0.7,0.7)-0.1)+0.015
for(i in seq(ppx))
{
  lines(px[1:3]+ppx[i], py[1:3]+ppy[i], col= "blue")
  lines(px[3:18]+ppx[i], py[3:18]+ppy[i], col= Cc[i])
  lines(px[18:33]+ppx[i], py[18:33]+ppy[i])
  lines(px[33:48]+ppx[i], py[33:48]+ppy[i], col= rev(Cc)[i])
  lines(px[48:51]+ppx[i], py[48:51]+ppy[i], col= "red")
  lines(px[51:100]+ppx[i], py[51:100]+ppy[i])
}

# RNA lines
t <- seq(0, 1, length=120)
p <- matrix(c(0,0, 0.15,2.25, 0.5,0.5, 0.85,-1.25, 1,1), nrow=5, ncol=2, byrow=TRUE)
p[,2] <-p[,2]/5
RNA <- bezier(t=t, p=p)*0.125
RNA <- RNA[11:110,]
x <- c(0.25, 0.27, 0.3, 0.35, 0.425, 0.46, 0.425)+0.09
y <- c(0.61,0.56,0.625,0.57,0.62,0.54,0.58)-0.01
for(i in seq(x))
{
  lines(RNA[1:5,1]+x[i], RNA[1:5,2]+y[i], col= "red")
  lines(RNA[5:35,1]+x[i], RNA[5:35,2]+y[i], col= Cc[i])
  lines(RNA[35:65,1]+x[i], RNA[35:65,2]+y[i])
  lines(RNA[65:95,1]+x[i], RNA[65:95,2]+y[i], col= rev(Cc)[i])
  lines(RNA[95:100,1]+x[i], RNA[95:100,2]+y[i], col= "blue")
}

# Sequencing
segments(0.79, 0.6, 0.81, 0.6, col= "red")
segments(0.81, 0.6, 0.86, 0.6, col= Cc[1])
segments(0.86, 0.6, 0.91, 0.6)
segments(0.91, 0.6, 0.96, 0.6, col= Cc[3])
segments(0.96, 0.6, 0.98, 0.6, col= "blue")
text(0.83, 0.6, "read 1", pos= 3, offset= 0.5, cex= 0.6)
text(0.94, 0.6, "read 2", pos= 3, offset= 0.5, cex= 0.6)
text(0.79, 0.55, "5'seq", pos= 4, offset= -0.5)
text(0.98, 0.55, "3'seq", pos= 2, offset= -0.5)

# Curved arrows enhancers to promoter
t <- seq(0, 1, length= 100)
p <- matrix(c(0.1765,0.865, 0.25,0.925, 0.365,0.86), nrow=3, ncol=2, byrow=TRUE)
lines(bezier(t=t, p=p))
p <- matrix(c(0.18,0.8675, 0.3,0.975, 0.55,0.86), nrow=3, ncol=2, byrow=TRUE)
lines(bezier(t=t, p=p))
polygon(c(0.17, 0.179, 0.187)-0.0025,
        c(0.86, 0.876, 0.867)-0.0025, border= NA, col= "black")

# Arrow plasmid to RNAs
segments(0.22,0.6,0.3,0.6)
polygon(c(0.3, 0.32, 0.3),
        c(0.6075, 0.6, 0.5925), border= NA, col= "black")

# Arrow RNAs to sequencing
segments(0.66,0.6,0.74,0.6)
polygon(c(0.74, 0.76, 0.74),
        c(0.6075, 0.6, 0.5925), border= NA, col= "black")

# Arrows sequencing
segments(0.79, 0.61, 0.82, 0.61)
polygon(c(0.82, 0.82, 0.835),
        c(0.6075, 0.615, 0.6075),
        col= "black",
        border= NA)
segments(0.95, 0.59, 0.98, 0.59)
polygon(c(0.95, 0.95, 0.935),
        c(0.5925, 0.585, 0.5925),
        col= "black",
        border= NA)

# 5' indidividual acitvities
pos <- seq(0, 0.07, length.out= 5)
for(i in seq(pos))
{
  segments(0.05, c(0.3-pos[i]), 0.12, c(0.3-pos[i]), col= adjustcolor("limegreen", 0.7), lwd= 3)
  segments(0.12, c(0.3-pos[i]), 0.19, c(0.3-pos[i]))
  segments(0.19, c(0.3-pos[i]), 0.26, c(0.3-pos[i]), col= grey.colors(8)[-1][i], lwd= 3)
}
rect(0.18, 0.32, 0.27, 0.21, lty= 3, lwd= 0.5)
text(0.225, 0.21, "control\nsequences", pos= 1, cex= 0.7, offset= 0.4)
text(0.155, 0.32, "5' individual\nactivity", pos= 3, offset= 0.5, col= adjustcolor("limegreen", 0.7))

# 3' indidividual acitvities
for(i in seq(pos))
{
  segments(0.4, c(0.3-pos[i]), 0.47, c(0.3-pos[i]), col= grey.colors(8)[-1][i], lwd= 3)
  segments(0.47, c(0.3-pos[i]), 0.54, c(0.3-pos[i]))
  segments(0.54, c(0.3-pos[i]), 0.61, c(0.3-pos[i]), col= adjustcolor("darkgoldenrod1", 0.7), lwd= 3)
}
rect(0.39, 0.32, 0.48, 0.21, lty= 3, lwd= 0.5)
text(0.435, 0.21, "control\nsequences", pos= 1, cex= 0.7, offset= 0.4)
text(0.505, 0.32, "3' individual\nactivity", pos= 3, offset= 0.5, col= adjustcolor("darkgoldenrod1", 0.7))

# Paired activity
segments(0.77, 0.3, 0.84, 0.3, col= adjustcolor("limegreen", 0.7), lwd= 3)
segments(0.84, 0.3, 0.91, 0.3)
segments(0.91, 0.3, 0.98, 0.3, col= adjustcolor("darkgoldenrod1", 0.7), lwd= 3)
text(0.875, 0.32, "5'/3' combined\nactitivy", pos= 3, offset= 0.5)

#------------------------------------------------#
# Luciferase
#------------------------------------------------#
# Backbone
segments(0, 0.05, 0.05, 0.05)
segments(0.3, 0.05, 0.7, 0.05)

# Promoter
segments(0.4, 0.05, 0.4, 0.1)
segments(0.4, 0.1, 0.45, 0.1, lend=2)
polygon(c(0.45,0.475,0.45),
        c(0.11,0.1,0.09), 
        border= NA, 
        col= "black")

# Boxes
rect(0.35, 0.025, 0.45, 0.075, col= "white")
text(0.4, 0.05, "CP")
rect(0.45, 0.025, 0.65, 0.075, col= "white")
text(0.55, 0.05, "luciferase")

# Enhancer pairs
Cc <- adjustcolor(vl_palette_few_categ(3), 0.7)
pos <- seq(0, 0.043, length.out= length(Cc))
for(i in seq(pos))
{
  segments(0.06, c(0.05-pos[i]), 0.13, c(0.05-pos[i]), col= Cc[i], lwd= 3)
  segments(0.13, c(0.05-pos[i]), 0.22, c(0.05-pos[i]))
  segments(0.22, c(0.05-pos[i]), 0.29, c(0.05-pos[i]), col= rev(Cc)[i], lwd= 3)
}

dev.off()