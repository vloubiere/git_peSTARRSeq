setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")

# Make matrix and clean
mat <- dcast(dat, -indL+L~indR+R, value.var = "log2FoldChange", sep = "__")
mat <- as.matrix(mat[, -1], 1)
colnames(mat) <- unlist(tstrsplit(colnames(mat), "__", keep= 2))
while(sum(is.na(mat))>0.05*nrow(mat)*ncol(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}
hm <- vl_heatmap(mat, 
                 cluster_rows= F,
                 cluster_cols= F, 
                 breaks = seq(-2, 8, 2),
                 col= viridis::viridis(6), 
                 show_rownames= F,
                 show_colnames= F, 
                 legend_title= "Activity (log2)")
Cc <- circlize::colorRamp2(hm$breaks, hm$col)

pdf("pdf/draft/heatmap_ordered_ind_act.pdf", width = 5, height = 5)
par(mar= c(5.1, 4.1, 4.1, 5.5),
    mgp= c(3, 0.5, 0))
plot(hm)
left <- matrix(dat[L %in% rownames(mat)][order(-indL), unique(indL)], ncol= 1)
rasterImage(Cc(left), 
            par("usr")[1]-strwidth("M")*2.5,
            par("usr")[3],
            par("usr")[1]-strwidth("M"),
            par("usr")[4], 
            xpd= NA)
right <- matrix(dat[R %in% colnames(mat)][order(indR), unique(indR)], nrow= 1)
rasterImage(Cc(right), 
            par("usr")[1],
            par("usr")[3]-strheight("M")*2.5,
            par("usr")[2],
            par("usr")[3]-strheight("M"), 
            xpd= NA)
title(xlab= "3' individual activity",
      ylab= "5' individual activity")
segments(par("usr")[1]-strwidth("M")*2.5,
         which.min(rev(abs(left-1))),
         par("usr")[2],
         which.min(rev(abs(left-1))),
         col= "white",
         xpd= NA)
segments(which.min(abs(right-1)),
         par("usr")[3]-strheight("M")*2.5,
         which.min(abs(right-1)),
         par("usr")[4],
         col= "white",
         xpd= NA)
text(par("usr")[1]-strwidth("M")*2.5,
     which.min(rev(abs(left-1))),
     labels = 1,
     xpd= NA, 
     pos= 2)
text(which.min(abs(right-1)),
     par("usr")[3]-strheight("M")*2.5,
     labels = 1,
     xpd= NA, 
     pos= 1)
dev.off()