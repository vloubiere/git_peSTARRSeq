setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("Rdata/vlli002_final_table.rds")
dat <- dat[!is.na(PC1L) & !is.na(PC1R)] # The PCA was only made for active pairs
dat[, L:= factor(L, unique(L[order(-PC1L)]))]
dat[, R:= factor(R, unique(R[order(PC1R)]))]

#---------------------------------------------#
# Order residuals matrix using PCA
#---------------------------------------------#
# Order and plot
im <- dcast(dat,
            L~R,
            value.var = "residuals")
hm <- vl_heatmap(as.matrix(im, 1), 
                 cluster_rows= F,
                 cluster_cols= F, 
                 breaks = c(-2, 0, 2),
                 show_rownames= F,
                 show_colnames= F, 
                 legend_title= "Residuals (log2)", 
                 show_legend = T,
                 plot = F)

#---------------------------------------------#
# Left side
#---------------------------------------------#
leftAct <- dat[, .(ind= indL[1]), keyby= L]
hm1 <- vl_heatmap(as.matrix(leftAct, 1),
                  cluster_rows= F,
                  cluster_cols= F, 
                  breaks = c(0, 7),
                  col = c("blue", "yellow"),
                  show_rownames= F,
                  show_colnames= F, 
                  legend_title= "Residuals (log2)", 
                  show_legend = F,
                  plot = F)
leftMot <- dat[, .SD[1, .(`AP-1_L`, Trl_L, twist_L, Dref_L)], keyby= L]
setnames(leftMot, function(x) gsub("_L", "", x))
hm2 <- vl_heatmap(as.matrix(leftMot, 1),
                  cluster_rows= F,
                  cluster_cols= F, 
                  breaks = c(0, 3),
                  col = c("white", "black"),
                  show_rownames= F,
                  show_colnames= T,
                  legend_title= "Residuals (log2)", 
                  show_legend = F,
                  plot = F)

#---------------------------------------------#
# Right side
#---------------------------------------------#
rightAct <- dat[, .(ind= indR[1]), keyby= R]
hm3 <- vl_heatmap(t(as.matrix(rightAct, 1)),
                  cluster_rows= F,
                  cluster_cols= F, 
                  breaks = c(0, 7),
                  col = c("blue", "yellow"),
                  show_rownames= F,
                  show_colnames= F, 
                  legend_title= "Residuals (log2)", 
                  show_legend = F,
                  plot = F)
rightMot <- dat[, .SD[1, .(Dref_R, twist_R, Trl_R, `AP-1_R`)], keyby= R]
setnames(rightMot, function(x) gsub("_R", "", x))
hm4 <- vl_heatmap(t(as.matrix(rightMot, 1)),
                  cluster_rows= F,
                  cluster_cols= F, 
                  breaks = c(0, 3),
                  col = c("white", "black"),
                  show_rownames= T,
                  show_colnames= F,
                  legend_title= "Residuals (log2)", 
                  show_legend = F,
                  plot = F)

#---------------------------------------------#
# Plot
#---------------------------------------------#
pdf("pdf/draft/heatmap_residuals_pca_order.pdf",
    width = 6,
    height = 5.25)
mat <- matrix(1:9, ncol= 3)
layout(mat,
       widths = c(0.175, 0.055, 1),
       heights = c(1, 0.055, 0.175))
par(mar= rep(.2, 4),
    oma= c(4,4,3,9),
    mgp= c(0.5,0.25,0),
    cex.lab= 1.5,
    las= 2)
plot(hm2)
title(ylab= "5' enhancer",
      xpd= NA)
plot.new()
plot.new()
plot(hm1)
plot.new()
plot.new()
plot(hm)
mar.lw <- diff(grconvertX(c(0,1), "line", "user"))
mar.lh <- diff(grconvertY(c(0,1), "line", "user"))
vl_heatkey(c(0,7), 
           col = c("blue", "yellow"),
           left = par("usr")[2]+mar.lw,
           top = par("usr")[4]-mar.lh*10,
           height= mar.lh*6,
           width= mar.lw, 
           main = "Individual act. (log2)")
vl_heatkey(0:3,
           col = colorRampPalette(c("white", "black"))(4),
           continuous = F,
           left = par("usr")[2]+mar.lw,
           top = par("usr")[4]-mar.lh*18,
           height= mar.lh*4,
           width= mar.lw, 
           main = "Motif counts")
plot(hm3)
par(las= 1)
plot(hm4)
title(xlab= "3' enhancer", 
      xpd= NA)
dev.off()
