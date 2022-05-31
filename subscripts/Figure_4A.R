setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(diagram)
require(bezier)
require(plotrix)
require(vlfunctions)

pdf("pdf/draft/Figure_4A.pdf", width = 4, height = 4)
par(mar= c(0,0,0,0),
    xaxt= "n",
    yaxt= "n",
    lend= 3,
    lwd= 1.3,
    cex= 1)

# Plot init
plot.new()

# Backbone
segments(0, 0.8, 0.25, 0.8)
segments(0.775, 0.8, 1, 0.8)

# Promoter
segments(0.075, 0.8, 0.075, 0.85)
segments(0.075, 0.85, 0.125, 0.85, lend=2)
polygon(c(0.125,0.15,0.125),
        c(0.86,0.85,0.84), 
        border= NA, 
        col= "black")

# Boxes
rect(0.025, 0.775, 0.125, 0.825, col= "white")
text(0.075, 0.8, "CP")
rect(0.825, 0.775, 0.975, 0.825, col= "white")
text(0.9, 0.8, "pA site")

# Enhancer pairs
Cc <- adjustcolor(vl_palette_few_categ(7), 0.7)
pos <- seq(0, 0.1, length.out= length(Cc))
for(i in seq(pos))
{
  segments(0.275, c(0.8-pos[i]), 0.35, c(0.8-pos[i]), col= Cc[i], lwd= 3)
  segments(0.375, c(0.8-pos[i]), 0.65, c(0.8-pos[i]), lty= "11")
  segments(0.675, c(0.8-pos[i]), 0.75, c(0.8-pos[i]), col= rev(Cc)[i], lwd= 3)
  segments(0.35, c(0.8-pos[i]), 0.375, c(0.8-pos[i]), col= "red")
  segments(0.65, c(0.8-pos[i]), 0.675, c(0.8-pos[i]), col= "blue")
}
text(0.3125, 0.805, "5'seq", pos= 3, offset= 0.25)
text(0.7125, 0.805, "3'seq", pos= 3, offset= 0.25)
text(0.5125, 0.805, "2kb spacer", pos= 3, offset= 0.25)

# Curved arrows enhancers to promoter
t <- seq(0, 1, length= 100)
p <- matrix(c(0.1765,0.865, 0.25,0.9, 0.30,0.86), nrow=3, ncol=2, byrow=TRUE)
lines(bezier(t=t, p=p))
p <- matrix(c(0.18,0.8675, 0.4,0.975, 0.7,0.86), nrow=3, ncol=2, byrow=TRUE)
lines(bezier(t=t, p=p))
polygon(c(0.17, 0.182, 0.187)-0.0025,
        c(0.862, 0.875, 0.866)-0.0025, border= NA, col= "black")

# Plasmids
px <- cos(seq(0, 2*pi, length.out= 100))
px <- px-min(px)
px <- px/max(px)/20
py <- sin(seq(0, 2*pi, length.out= 100))
py <- py-min(py)
py <- py/max(py)/20
set.seed(1)
ppx <- jitter(c(0.1,0.175, 0.1375, 0.2125, 0.1,0.175))-0.1
set.seed(8)
ppy <- jitter(c(0.6,0.6,0.65,0.65,0.7,0.7,0.7)-0.1)+0.015
for(i in seq(ppx))
{
  lines(px[1:12]+ppx[i], py[1:12]+ppy[i], col= Cc[i])
  lines(px[12:15]+ppx[i], py[12:15]+ppy[i], col= "blue")
  lines(px[15:60]+ppx[i], py[15:60]+ppy[i], lty= "11")
  lines(px[60:63]+ppx[i], py[60:63]+ppy[i], col= "red")
  lines(px[63:75]+ppx[i], py[63:75]+ppy[i], col= rev(Cc)[i])
  lines(px[75:100]+ppx[i], py[75:100]+ppy[i])
}

# RNA lines
t <- seq(0, 1, length=120)
p <- matrix(c(0,0, 0.15,2.25, 0.5,0.5, 0.85,-1.25, 1,1), nrow=5, ncol=2, byrow=TRUE)
p[,2] <-p[,2]/5
RNA <- bezier(t=t, p=p)*0.125
RNA <- RNA[11:110,]
x <- c(0.335, 0.35, 0.425, 0.375, 0.40)-0.1
y <- c(0.625,0.57,0.62,0.54,0.58)-0.01
for(i in seq(x))
{
  lines(RNA[1:20,1]+x[i], RNA[1:20,2]+y[i], col= Cc[i])
  lines(RNA[20:25,1]+x[i], RNA[20:25,2]+y[i], col= "red")
  lines(RNA[25:65,1]+x[i], RNA[25:65,2]+y[i], lty= "11")
  lines(RNA[65:70,1]+x[i], RNA[65:70,2]+y[i], col= "blue")
  lines(RNA[70:95,1]+x[i], RNA[70:95,2]+y[i], col= rev(Cc)[i])
}

# Reverse PCR
px <- cos(seq(0, 2*pi, length.out= 100))
px <- px-min(px)
px <- px/max(px)/30
py <- sin(seq(0, 2*pi, length.out= 100))
py <- py-min(py)
py <- py/max(py)/30
set.seed(1)
ppx <- jitter(c(0.1,0.15, 0.13, 0.18))+0.43
set.seed(8)
ppy <- jitter(c(0.6,0.6,0.65,0.65,0.7)-0.1)+0.06
for(i in seq(ppx))
{
  lines(px[1:5]+ppx[i], py[1:5]+ppy[i], col= "blue")
  lines(px[5:15]+ppx[i], py[5:15]+ppy[i], col= Cc[i])
  lines(px[15:35]+ppx[i], py[15:35]+ppy[i], col= "black")
  lines(px[35:45]+ppx[i], py[35:45]+ppy[i], col= rev(Cc)[i])
  lines(px[45:49]+ppx[i], py[45:49]+ppy[i], col= "red")
  lines(px[50:100]+ppx[i], py[50:100]+ppy[i], lty= "11")
  segments(px[75]-0.005+ppx[i],
           py[75]-0.0025+ppy[i],
           px[75]+0.0035+ppx[i],
           py[75]-0.0175+ppy[i],
           lwd= 0.5)
  segments(px[75]+0.005+ppx[i],
           py[75]-0.0025+ppy[i],
           px[75]-0.0035+ppx[i],
           py[75]-0.0175+ppy[i],
           lwd= 0.5)
  points(px[75]-0.005+ppx[i],
         py[75]-0.02+ppy[i], 
         cex= 0.3,
         lwd= 0.5)
  points(px[75]+0.005+ppx[i],
         py[75]-0.02+ppy[i], 
         cex= 0.3,
         lwd= 0.5)
}


# Sequencing
segments(0.77, 0.6, 0.795, 0.6, lty= "11")
segments(0.795, 0.6, 0.81, 0.6, col= "red")
segments(0.81, 0.6, 0.86, 0.6, col= Cc[1])
segments(0.86, 0.6, 0.91, 0.6)
segments(0.91, 0.6, 0.96, 0.6, col= Cc[3])
segments(0.96, 0.6, 0.975, 0.6, col= "blue")
segments(0.98, 0.6, 1, 0.6, lty= "11")
text(0.83, 0.6, "read 1", pos= 3, offset= 0.5, cex= 0.6)
text(0.94, 0.6, "read 2", pos= 3, offset= 0.5, cex= 0.6)
text(0.79, 0.55, "5'seq", pos= 4, offset= -0.5)
text(0.98, 0.55, "3'seq", pos= 2, offset= -0.5)
segments(0.795, 0.61, 0.825, 0.61)
polygon(c(0.825, 0.825, 0.84),
        c(0.6075, 0.615, 0.6075),
        col= "black",
        border= NA)
segments(0.945, 0.59, 0.975, 0.59)
polygon(c(0.945, 0.945, 0.93),
        c(0.5925, 0.585, 0.5925),
        col= "black",
        border= NA)


# Arrow plasmid to RNAs
segments(0.19,0.6,0.22,0.6)
polygon(c(0.22, 0.23, 0.22),
        c(0.605, 0.6, 0.595), 
        border= NA,
        col= "black")

# Arrow RNAs to ligation
segments(0.46,0.6,0.49,0.6)
polygon(c(0.49, 0.5, 0.49),
        c(0.605, 0.6, 0.595), 
        border= NA, 
        col= "black")
text(0.475, 
     0.6, 
     "Ligation", 
     pos= 3,
     offset= 0.5,
     cex= 0.6)

# Arrow ligation to sequencing
segments(0.68,0.6,0.71,0.6)
polygon(c(0.71, 0.72, 0.71),
        c(0.605, 0.6, 0.595), 
        border= NA, 
        col= "black")


dev.off()