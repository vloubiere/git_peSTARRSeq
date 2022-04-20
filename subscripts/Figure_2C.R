setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
cl <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")
dat <- melt(as.data.table(cl$x, keep.rownames = "L"), id.vars= "L", variable.name = "R")
dat[cl$rows, cl_L:= paste0("5' ", i.cl), on= "L==name"]
dat[cl$cols, cl_R:= paste0("3' ", i.cl), on= "R==name"]
dat[, cl:= paste0(cl_L, "\n", cl_R)]
Cc <- circlize::colorRamp2(seq(-3,3, length.out= 20),
                           vl_palette_blueWhiteRed(20))

pdf("pdf/draft/Figure_2C.pdf", 
    width = 3,
    height = 4.5)
par(las= 2)
vl_boxplot(value~cl, 
           dat, 
           ylab= "Observed/Additive (log2)",
           violin= T, 
           violcol= adjustcolor(Cc(dat[, median(value, na.rm= T), keyby= cl]$V1), 0.7))
abline(h=0, 
       lty= 2)
dev.off()