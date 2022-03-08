setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(vioplot)
require(bezier)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- lib[vllib=="vllib002"]
dat <- feat$add_feature(dat, feat$lib)
dat <- feat$add_feature(dat, feat$closest_TSSs)

enh <- unique(dat[, c(L, R)])
pairs <- data.table(A= enh)
pairs <- pairs[, .(B= enh[-(1:.GRP)]), A]
pairs[dat, c("log2FoldChange.x", "group.x"):= .(i.log2FoldChange, i.group_L), on= c("A==L", "B==R")]
pairs[dat, c("log2FoldChange.y", "group.y"):= .(i.log2FoldChange, i.group_L), on= c("A==R", "B==L")]
pairs <- na.omit(pairs)
pairs[, col:= fcase(group.x=="control" & group.y=="control", "black",
                    group.x!="control" & group.y!="control", "#74C27A",
                    default= "royalblue2")]
set.seed(1)
pairs <- pairs[sample(nrow(pairs), nrow(pairs))]

pdf("pdf/draft/Figure_1D.pdf", 4.5, 4.5)
par(mar= c(5,5,3,3),
    las= 1)
plot(pairs$log2FoldChange.x,
     pairs$log2FoldChange.y,
     col= adjustcolor(pairs$col, 0.5),
     xlab= "Activity E1/E2 (log2)",
     ylab= "Activity E2/E1 (log2)",
     cex= 0.5,
     pch= 16,
     xlim= c(-4.2, 10.1),
     ylim= c(-4.2, 10.1))

# Add curved arrow
px <- cos(seq(pi, pi*3.5, length.out= 250))
px <- px-min(px)
px <- px-8.5
py <- sin(seq(pi, pi*3.5, length.out= 250))
py <- py-min(py)
py <- py-8.5
t <- seq(0, 1, length= 250)
p <- cbind(px, py)
p[51:250,1] <- p[51:250,1]+0.2
p[1:200,2] <- p[1:200,2]+0.2
p <- rbind(c(-8.5,-2),
           p,
           c(-2,-8.5))
lines(bezier(t=t, p=p), 
      lwd= 2, 
      xpd= T)
polygon(c(-8.3, -8.5, -8.7),
        c(-2, -1.6, -2),
        border= NA,
        col= "black",
        xpd= T)
polygon(c(-2, -1.6, -2),
        c(-8.3, -8.5, -8.7),
        border= NA,
        col= "black",
        xpd= T)

# Pair cartoons
segments(c(-6,-4.5,-3), 
         -7.75, 
         c(-4.5,-3,-1.5), 
         -7.75, 
         col= c("#E69F00", "black", "#CD748B"), 
         lwd= c(6,3,6),
         lend= 3,
         xpd= T)
segments(-7.75, 
         c(-6,-4.5,-3),
         -7.75,
         c(-4.5,-3,-1.5), 
         col= c("#CD748B", "black", "#E69F00"), 
         lwd= c(6,3,6),
         lend= 3,
         xpd= T)

# Densities
bottom <- par("usr")[4]
top <- grconvertY(1, "nfc", "user")-strheight("M")/2
x_dens <- pairs[, density(log2FoldChange.x)[c("x", "y")], col]
x_dens[, adj.y:= bottom+y/max(y)*(top-bottom), col]
x_dens[, polygon(x, adj.y, xpd= T, border=NA, col= adjustcolor(col[1], 0.5)), col]

left <- par("usr")[2]
right <- grconvertX(1, "nfc", "user")-strwidth("M")/2
y_dens <- pairs[, setNames(density(log2FoldChange.y)[c("x", "y")], c("y", "x")), col]
y_dens[, adj.x:= left+x/max(x)*(right-left), col]
y_dens[, polygon(adj.x, y, xpd= T, border=NA, col= adjustcolor(col[1], 0.5)), col]

legend("topleft", 
       col= adjustcolor(c("royalblue2","#74C27A","black"), 0.8),
       legend= c(paste0("enh./ctl. OR ctl./enh. (PCC=", 
                        round(cor.test(pairs[col=="royalblue2", log2FoldChange.x],
                                       pairs[col=="royalblue2", log2FoldChange.y])$estimate, 2), ")"), 
                 paste0("enh./enh. (PCC=", 
                        round(cor.test(pairs[col=="#74C27A", log2FoldChange.x],
                                       pairs[col=="#74C27A", log2FoldChange.y])$estimate, 2), ")"),
                 "ctl./ctl."),
                 bty= "n",
                 pch= 16,
       cex= 0.7)
abline(0, 1, lty=2)
dev.off()