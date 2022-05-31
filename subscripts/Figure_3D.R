setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
cl <- readRDS("Rdata/vllib016_clustering_additive_scores_draft_figure.rds")
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib016" 
                                                & class== "enh./enh."
                                                & grepl("dev|hk", L)
                                                & grepl("dev|hk", R)]
dat[, diff:= log2FoldChange-additive]
dat[cl$rows, cl_L:= i.cl, on= "L==name"]
dat[cl$cols, cl_R:= i.cl, on= "R==name"]
dat[, cl:= paste0("5'", cl_L, "-3'", cl_R)]
dat[, ord:= median(diff), cl]
setorderv(dat, "ord", -1)
dat[, cl:= factor(cl, unique(dat[, cl]))]
Cc <- circlize::colorRamp2(cl$breaks,
                           cl$col)

#----------------------------#
# PLOT
#----------------------------#
pdf("pdf/draft/Figure_3D.pdf",
    height= 2.3,
    width= 2)
par(las= 2,
    mar= c(3.5,3,0.5,1),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2)
vl_boxplot(split(dat$diff, dat$cl),
           violin= T,
           violcol= adjustcolor(Cc(sapply(split(dat$diff, dat$cl), median, na.rm= T)), 0.7),
           ylab= "Obs./Exp. Add. (log2)", 
           tilt.names = T)
abline(h=0, 
       lty= 2)
dev.off()