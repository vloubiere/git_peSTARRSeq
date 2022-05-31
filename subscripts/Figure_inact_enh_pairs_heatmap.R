setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class=="inact./enh."]
feat <- readRDS("Rdata/uniq_enh_feat/lib_genomic_dat.rds")
dat[, diff:= log2FoldChange-additive]

mat <- as.matrix(dcast(dat, L~R, value.var= "diff"), 1)
mat <- mat[apply(mat, 1, function(x) sum(is.na(x))<=250),apply(mat, 2, function(x) sum(is.na(x))<=250)]
cl <- vl_heatmap(mat,
                 cutree_rows = 2, 
                 cutree_cols = 2,
                 breaks = seq(-3, 3, length.out= 30), 
                 col = vl_palette_blueWhiteRed(30, rep_white = 5),
                 legend_title = "Obs./Add. (log2)", 
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
cl$rows[, cl:= as.character(cl)]
cl$rows[cl==1, c("cl", "col"):= .("A", "grey90")]
cl$rows[cl==2, c("cl", "col"):= .("B", "grey60")]
cl$cols[, cl:= as.character(cl)]
cl$cols[cl==1, c("cl", "col"):= .("A", "grey90")]
cl$cols[cl==2, c("cl", "col"):= .("B", "grey60")]
cl$rows[unique(dat[,.(name=L, median_L)]), ind_act:= median_L, on= "name"]
ind_L <- cl$rows[(order), ind_act]
cl$cols[unique(dat[,.(name=R, median_R)]), ind_act:= median_R, on= "name"]
ind_R <- cl$cols[(order), ind_act]
# Add motif enrichment
cl$mot_enr_L$seq <- lib[cl$rows$name, oligo_full_sequence, on= "ID_vl"]
cl$mot_enr_L$counts <- vl_motif_counts(cl$mot_enr_L$seq)
cl$mot_enr_L$counts <- cl$mot_enr_L$counts[, apply(cl$mot_enr_L$counts, 2, function(x) sum(x>0))>20]
cl$mot_enr_L$enr <- vl_motif_enrich(cl$mot_enr_L$counts[cl$rows$cl=="A",],
                                    cl$mot_enr_L$counts[cl$rows$cl=="B",],
                                    plot= F)
cl$mot_enr_R$seq <- lib[cl$cols$name, oligo_full_sequence, on= "ID_vl"]
cl$mot_enr_R$counts <- vl_motif_counts(cl$mot_enr_R$seq)
cl$mot_enr_R$counts <- cl$mot_enr_R$counts[, apply(cl$mot_enr_R$counts, 2, function(x) sum(x>0))>20]
cl$mot_enr_R$enr <- vl_motif_enrich(cl$mot_enr_R$counts[cl$cols$cl=="A",],
                                    cl$mot_enr_R$counts[cl$cols$cl=="B",],
                                    plot= F)
saveRDS(cl, "Rdata/vllib002_clustering_inactive_pairs_additive_scores_draft_figure.rds")

pdf("pdf/draft/Figure_inact_enh_pairs_heatmap.pdf", 
    width= 8, 
    height = 7)
par(mar= c(3.5,3.5,3,10),
    mgp= c(3,0.15,0),
    tcl= -0.2)

# Heatmap
plot(cl)

# Left individual activities
right <- par("usr")[1]-strwidth("M", cex= 0.5)
width <- right-grconvertX(1, "line", "user")
rect(right-(ind_L/max(ind_L)*width), 
     rev(seq(ind_L))+0.5, 
     right, 
     rev(seq(ind_L))-0.5, 
     border= NA, 
     xpd= T, 
     col= adjustcolor("#0C3A0E", 0.7))
ticks <- axisTicks(c(0, max(ind_L)), log= F)
at <- right-(ticks/max(ind_L)*width)
axis(3, 
     at = range(at), 
     labels = range(ticks),
     xpd= T, 
     line= 0.25, 
     cex.axis= 0.5)
text(mean(at),
     par("usr")[4],
     "5' Individual\nact. (log2)",
     xpd= T, 
     pos= 3, 
     offset= 1.25,
     cex= 0.7)
text(grconvertX(0.5, "line", "user"),
     mean(par("usr")[c(3,4)]),
     "5' enhancer",
     srt= 90,
     xpd= T)

# Right individual activities
par(mgp= c(3,0.35,0))
bottom <- grconvertY(1, "line", "user")
height <- par("usr")[3]-strheight("M")*0.5-bottom
rect(seq(ind_R)+0.5,
     bottom+(ind_R/max(ind_R)*height), 
     seq(ind_R)-0.5,
     bottom,
     border= NA, 
     xpd= T, 
     col= adjustcolor("#0C3A0E", 0.7))
ticks <- axisTicks(c(0, max(ind_R)), log= F)
at <- bottom+(ticks/max(ind_R)*height)
axis(2, 
     at = range(at), 
     labels = range(ticks),
     xpd= T, 
     line= 0.25, 
     cex.axis= 0.5,
     las= 2)
text(par("usr")[1],
     mean(at),
     "3' Individual\nact. (log2)",
     xpd= T,
     pos= 2, 
     offset= 0.5,
     cex= 0.7)
text(mean(par("usr")[c(1,2)]),
     grconvertY(0.5, "line", "user"),
     "3' enhancer",
     xpd= T)
dev.off()