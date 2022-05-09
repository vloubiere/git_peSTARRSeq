setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib016" 
                                                & class== "enh./enh."
                                                & grepl("dev|hk", L)
                                                & grepl("dev|hk", R)]
feat <- readRDS("Rdata/uniq_enh_feat/lib_genomic_dat.rds")
dat[, diff:= log2FoldChange-additive]

mat <- as.matrix(dcast(dat, L~R, value.var= "diff"), 1)
cl <- vl_heatmap(mat, 
                 cutree_rows = 3, 
                 cutree_cols = 2,
                 breaks = seq(-3, 3, length.out= 30), 
                 col = vl_palette_blueWhiteRed(30, rep_white = 5),
                 legend_title = "Obs./Add. (log2)",
                 clustering_method = "ward.D",
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
cl$rows[, cl:= as.character(cl)]
cl$rows[cl==1, c("cl", "col"):= .("A", "grey90")]
cl$rows[cl==3, c("cl", "col"):= .("B", "grey65")]
cl$rows[cl==2, c("cl", "col"):= .("C", "grey40")]
cl$rows[unique(dat[,.(name=L, median_L)]), ind_act:= median_L, on= "name"]
ind_L <- cl$rows[(order), ind_act]
cl$cols[, cl:= as.character(cl)]
cl$cols[cl==1, c("cl", "col"):= .("A", "grey90")]
cl$cols[cl==2, c("cl", "col"):= .("B", "grey60")]
cl$cols[unique(dat[,.(name=R, median_R)]), ind_act:= median_R, on= "name"]
ind_R <- cl$cols[(order), ind_act]
# Add motif enrichment
cl$rows[lib, seq:= i.oligo_full_sequence, on= "name==ID"]
cl$cols[lib, seq:= i.oligo_full_sequence, on= "name==ID"]
sequences <- data.table(rbind(cl$rows[, .(cl= paste0("5' ", cl), seq)],
                              cl$cols[, .(cl= paste0("3' ", cl), seq)],
                              lib[grepl("^control", ID), .(cl= "control", seq= oligo_full_sequence)]))
cl$mot$counts <- vl_motif_counts(sequences$seq)
cl$mot$counts <- cl$mot$counts[, apply(cl$mot$counts, 2, function(x) sum(x, na.rm= T)>10)]
cl$mot$enr <- vl_motif_cl_enrich(cl$mot$counts,
                                 cl_IDs = sequences$cl, 
                                 control_cl = "control",
                                 plot= F)
saveRDS(cl,
        "Rdata/vllib016_clustering_additive_scores_draft_figure.rds")

pdf("pdf/draft/Figure_3C.pdf", 
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
Cc <- fcase(grepl("^dev", cl$rows$name), "#74C27A", 
            grepl("^hk", cl$rows$name), "tomato")
Cc <- rev(Cc[cl$rows$order])
rect(xleft = par("usr")[1]-strwidth("M")*2,
     ybottom = seq(nrow(cl$rows)-1)-0.5,
     xright = par("usr")[1]-strwidth("M")*0.5,
     ytop = seq(nrow(cl$rows))[-1]+0.5,
     xpd= T,
     border= NA,
     col= Cc)
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
Cc <- fcase(grepl("^dev", cl$cols$name), "#74C27A", 
            grepl("^hk", cl$cols$name), "tomato")
Cc <- Cc[cl$cols$order]
rect(xleft = seq(nrow(cl$cols)-1)-0.5,
     ybottom = par("usr")[3]-strheight("M")*2,
     xright = seq(nrow(cl$cols))[-1]+0.5,
     ytop = par("usr")[3]-strheight("M")*0.5,
     xpd= T,
     border= NA,
     col= Cc)
text(mean(par("usr")[c(1,2)]),
     grconvertY(0.5, "line", "user"),
     "3' enhancer",
     xpd= T)
legend(par("usr")[2], 
       par("usr")[3]-strheight("M")*0.5, 
       fill= c("#74C27A", "tomato"), 
       legend= c("Developmental", "Housekeeping"), 
       xpd= T,
       bty= "n")
dev.off()
