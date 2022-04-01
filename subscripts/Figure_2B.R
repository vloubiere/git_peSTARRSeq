setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class== "enh./enh."]
feat <- readRDS("Rdata/uniq_enh_feat/lib_genomic_dat.rds")
dat[, diff:= log2FoldChange-additive]

mat <- as.matrix(dcast(dat, L~R, value.var= "diff"), 1)
cl <- vl_heatmap(mat, 
                 cutree_rows = 2, 
                 cutree_cols = 2,
                 breaks = seq(-3,3, length.out= 20), 
                 col = vl_palette_blueWhiteRed(20),
                 legend_title = "Obs./Add. (log2)", 
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
cl$rows[unique(dat[,.(name=L, median_L)]), ind_act:= median_L, on= "name"]
ind_L <- cl$rows[(order), ind_act]
cl$cols[unique(dat[,.(name=R, median_R)]), ind_act:= median_R, on= "name"]
ind_R <- cl$cols[(order), ind_act]
# Add motif enrichment
cl$mot_enr_L$seq <- lib[cl$rows$name, oligo_full_sequence, on= "ID_vl"]
cl$mot_enr_L$counts <- vl_motif_counts(cl$mot_enr_L$seq)
cl$mot_enr_L$enr <- vl_motif_enrich(cl$mot_enr_L$counts[cl$rows$cl==2,],
                                    cl$mot_enr_L$counts[cl$rows$cl==1,],
                                    plot= F)
cl$mot_enr_R$seq <- lib[cl$cols$name, oligo_full_sequence, on= "ID_vl"]
cl$mot_enr_R$counts <- vl_motif_counts(cl$mot_enr_R$seq)
cl$mot_enr_R$enr <- vl_motif_enrich(cl$mot_enr_R$counts[cl$cols$cl==2,],
                                    cl$mot_enr_R$counts[cl$cols$cl==1,],
                                    plot= F)
saveRDS(cl, "Rdata/vllib002_clustering_additive_scores_draft_figure.rds")

pdf("pdf/draft/Figure_2B.pdf", 
    width= 10, 
    height = 9)
par(mar= c(6.5,6.5,3,10),
    mgp= c(3,0.15,0))

# Heatmap
plot(cl)

# Left individual activities
right <- par("usr")[1]-strwidth("M")*2.5
width <- right-grconvertX(1.5, "line", "user")
rect(right-(ind_L/max(ind_L)*width), 
     rev(seq(ind_L))+0.5, 
     right, 
     rev(seq(ind_L))-0.5, 
     border= NA, 
     xpd= T, 
     col= adjustcolor("#0C3A0E", 0.7))
ticks <- axisTicks(c(0,max(ind_L)), log= F)
at <- right-(ticks/max(ind_L)*width)
axis(3, 
     at = at, 
     labels = ticks,
     xpd= T, 
     line= 0.25, 
     cex.axis= 0.5,
     tck= -0.005)
text(mean(at),
     par("usr")[4],
     "5' Individual\nact. (log2)",
     xpd= T, 
     pos= 3, 
     offset= 1.5,
     cex= 0.7)
clr_pos <- cumsum(rev(table(cl$rows$cl)))
rect(xleft = par("usr")[1]-strwidth("M")*2,
     ybottom = c(1, clr_pos[-length(clr_pos)]),
     xright = par("usr")[1]-strwidth("M")*0.5,
     ytop = clr_pos,
     xpd= T,
     border= NA,
     col= adjustcolor(c("grey20", "grey70"), 0.7))
text(x = par("usr")[1]-strwidth("M")*1.25,
     y = clr_pos-diff(c(1, clr_pos))/2,
     c("5' cluster 2", "5' cluster 1"),
     xpd= T,
     srt= 90,
     col= c("white", "black"))
text(grconvertX(0.5, "line", "user"),
     mean(par("usr")[c(3,4)]),
     "5' enhancer",
     srt= 90,
     xpd= T)

# Right individual activities
par(mgp= c(3,0.35,0))
bottom <- grconvertY(1.5, "line", "user")
height <- par("usr")[3]-strheight("M")*2.5-bottom
rect(seq(ind_R)+0.5,
     bottom+(ind_R/max(ind_R)*height), 
     seq(ind_R)-0.5,
     bottom,
     border= NA, 
     xpd= T, 
     col= adjustcolor("#0C3A0E", 0.7))
ticks <- axisTicks(c(0,max(ind_R)), log= F)
at <- bottom+(ticks/max(ind_R)*height)
axis(2, 
     at = at, 
     labels = ticks,
     xpd= T, 
     line= 0.25, 
     cex.axis= 0.5,
     las= 2,
     tck= -0.005)
text(par("usr")[1],
     mean(at),
     "3' Individual\nact. (log2)",
     xpd= T,
     pos= 2, 
     offset= 1.25,
     cex= 0.7)
clc_pos <- cumsum(rev(table(cl$cols$cl)))
rect(xleft = c(1, clc_pos[-length(clc_pos)]),
     ybottom = par("usr")[3]-strheight("M")*2,
     xright = clc_pos,
     ytop = par("usr")[3]-strheight("M")*0.5,
     xpd= T,
     border= NA,
     col= adjustcolor(c("grey20", "grey70"), 0.7))
text(x = clc_pos-diff(c(1, clc_pos))/2,
     y = par("usr")[3]-strheight("M")*1.25,
     c("3' cluster 2", "3' cluster 1"),
     xpd= T,
     col= c("white", "black"))
text(mean(par("usr")[c(1,2)]),
     grconvertY(0.5, "line", "user"),
     "3' enhancer",
     xpd= T)
dev.off()