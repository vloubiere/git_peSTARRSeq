setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(vioplot)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")

# dat
dat <- lib[vllib=="vllib015"][act_wilcox_L<0.05 & act_wilcox_R<0.05]
dat <- feat$add_feature(dat, feat$lib)
dat <- dat[!group_L %in% c("control", "repressor") & !group_R %in% c("control", "repressor")]
dat <- feat$add_feature(dat, feat$top_motifs)
dat[, diff:= log2FoldChange-additive]
dat[, col:= colorRampPalette(unlist(.BY))(3)[2], .(col_L, col_R)]

#-----------------------------------------------#
# Computations
#-----------------------------------------------#
# Residuals clustering
obj <- vl_heatmap(as.matrix(dcast(dat, L~R, value.var = "diff"), 1), 
                  cutree_rows = 4, 
                  cutree_cols = 5, 
                  plot = F, 
                  clustering_method = "ward.D2")
dat[obj$result_DT[, .(rn, rn_cl)], cl_L:= i.rn_cl, on= "L==rn"]
dat[obj$result_DT[, .(variable, var_cl)], cl_R:= i.var_cl, on= "R==variable"]

# motif enrichment residual clusters
motL <- dat[, .SD[1], .(cl_L, L), .SDcols= patterns("top_motif.*_L")]
enrL <- vl_motif_cl_enrich(as.matrix(motL[, -1], 1), 
                           motL$cl_L, 
                           plot = F)
motR <- dat[, .SD[1], .(cl_R, R), .SDcols= patterns("top_motif.*_R")]
enrR <- vl_motif_cl_enrich(as.matrix(motR[, -1], 1), 
                           motR$cl_R, 
                           plot = F)

# Split between Jra Dref motifs
motif_split_func <- function(dat, 
                             prefixMot1,
                             prefixMot2,
                             name1,
                             name2)
{
  dat <- copy(dat)
  L1 <- dat[,get(paste0(prefixMot1, "_L"))]
  L2 <- dat[,get(paste0(prefixMot2, "_L"))]
  R1 <- dat[,get(paste0(prefixMot1, "_R"))]
  R2 <- dat[,get(paste0(prefixMot2, "_R"))]
  dat[, class_L:= fcase(L1>0 & L2==0, name1,
                        L1==0 & L2>0, name2,
                        L1>0 & L2>0, "both",
                        default= "none")]
  dat[, class_R:= fcase(R1>0 & R2==0, name1,
                        R1==0 & R2>0, name2,
                        R1>0 & R2>0, "both",
                        default= "none")]
  dat[, class_L:= factor(class_L, levels= c("both", name1, name2, "none"))]
  dat[, class_R:= factor(class_R, levels= c("none", name2, name1, "both"))]
  return(list(act= dcast(dat, class_L~class_R, value.var = "log2FoldChange", fun.aggregate = median),
              diff= dcast(dat, class_L~class_R, value.var = "diff", fun.aggregate = median)))
}

# Rank classes for residuals violin plots
dat[group_L %in% c("dev", "hk")
    & group_R %in% c("dev", "hk"), 
    c("res_group", "median_diff"):= .(paste0(group_L, "+", group_R), median(diff)), 
    .(group_L, group_R)]
dat[group_L=="dev" & group_R=="hk", col:= "royalblue"]
dat[group_L=="hk" & group_R=="dev", col:= "gold"]
setorderv(dat, "median_diff", -1)
dat[, res_group:= factor(res_group, levels = unique(res_group))]

#-----------------------------------------------#
# PLOT
#-----------------------------------------------#
Cc <- colorRampPalette(c("blue", "red"))(4)

pdf("pdf/analyses/DSCP_screen.pdf", 
    width = 18, 
    height = 15)
par(las= 1, 
    mfrow= c(3,3),
    mar= c(5.1, 18, 4.1, 18))

# Violin plot observed and expected
box <- vl_boxplot(dat[, .(Observed= log2FoldChange, Additive= additive)],
                  ylab= "Activity",
                  lty= 1, 
                  outline= T, 
                  violin= T, 
                  compute_pval= c(1,2))
abline(h= box$obj[1, Q2], lty= 2)

# Violin plot residuals
box <- vl_boxplot(dat[, diff],
                  ylab= "Observed/Additive (log2)",
                  lty= 1, 
                  outline= T,
                  violin= T, 
                  names= F)
abline(h= 0, lty= 2)

# quantiles
pol <- box[1, .(x= unlist(x), at, y= unlist(y))]
low_cutoff <- box[1,Q1]
high_cutoff <- box[1,Q3]
low <- pol[y<low_cutoff, .(x= c(x, x[.N]), y= c(y, low_cutoff), at= at[1])][, .(c(at+x, at-rev(x)), c(y, rev(y)))]
high <- pol[y>high_cutoff, .(x= c(x[1], x), y= c(high_cutoff, y), at= at[1])][, .(c(at+x, at-rev(x)), c(y, rev(y)))]
polygon(low, 
        border= NA, 
        col= adjustcolor("cornflowerblue", 0.5))
polygon(high, 
        border= NA, 
        col= adjustcolor("tomato", 0.5))

# Residuals clustering heatmap
par(mar= c(2, 2, 4.1, 12))
plot(obj, 
     auto_margins= F,
     cutree_rows= obj$cutree_rows,
     cutree_cols= obj$cutree_cols, 
     show_rownames = F,
     show_colnames = F, 
     breaks = c(-4,0,4), 
     legend_title= "Obs./Add. (log2)")

# Motif enrichment left enhancer
par(cex= 0.8,
    mar= c(2, 20, 2, 10))
vl_motif_cl_enrich_plot_only(obj = enrL,
                             x_breaks = c(1,2,3),
                             col = Cc,
                             padj_cutoff = 0.05,
                             N_top = 5,
                             auto_margins = F)

# Motif enrichment right enhancer
vl_motif_cl_enrich_plot_only(obj = enrR,
                             x_breaks = c(1,2,3),
                             col = Cc,
                             padj_cutoff = 0.05,
                             N_top = 5,
                             auto_margins = F)

# Balloon plot motifs
par(mar= c(5, 5, 2, 10))
cmb <- data.table(prefixMot1= c("top_motif_kay_Jra_1", "top_motif_Mad_4", "top_motif_grn_3"),
                  prefixMot2= c("top_motif_Dref", "top_motif_Dref", "top_motif_HLH3B"),
                  name1= c("kay/Jra", "Mad", "grn"),
                  name2= c("Dref", "Dref", "HLH3B"))
cmb[, {
  .c <- motif_split_func(dat, 
                         prefixMot1,
                         prefixMot2, 
                         name1,
                         name2)
  vl_balloons_plot(x = as.matrix(.c$act, 1),
                   x_breaks = c(2,3,4,5),
                   color_var = as.matrix(.c$diff, 1), 
                   color_breaks = c(-1,-0.2,0.2,1),
                   balloon_size_legend = "Activity (log2)",
                   balloon_col_legend = "Obs./Add. (log2)",
                   cex.balloons = 1.3,
                   auto_margins = F,
                   col = c("cornflowerblue", "white", "white", "tomato"))
  print("")
}, (cmb)]

# Additive model
par(mar= c(8, 8, 2, 8))
smoothScatter(dat[, .(additive, log2FoldChange)],
              colramp = colorRampPalette(c("white", "darkgrey")),
              ylab= "Activity (log2)",
              xlab= "Additive (log2)",
              xlim= c(1, 8),
              ylim= c(-1, 9))
dat[SJ(group_L= c("dev", "hk", "dev", "hk"),
       group_R= c("dev", "hk", "hk", "dev")), {
      points(additive, 
             log2FoldChange,
             col= adjustcolor(col[1], 0.5),
             pch= 19,
             cex= 0.5)
         if(group_R=="hk")
           abline(lm(log2FoldChange~additive), 
                  col= col[1],
                  lty= 2)
      
    }, by= .(group_L, group_R),
    on= c("group_L", "group_R")]
abline(0,1)

# Additive model MA plot style
smoothScatter(dat[, .(additive, log2FoldChange-additive)],
              colramp = colorRampPalette(c("white", "darkgrey")),
              xlab= "Additive (log2)",
              ylab= "Obs./Add. (log2)")
dat[SJ(group_L= c("dev", "hk", "dev", "hk"),
       group_R= c("dev", "hk", "hk", "dev")), {
         points(additive, 
                log2FoldChange-additive,
                col= adjustcolor(col[1], 0.5),
                pch= 19,
                cex= 0.25)
       }, by= .(group_L, group_R),
    on= c("group_L", "group_R")]
abline(h= 0)

# Violin plot residuals
yl <- par("usr")[3:4]
par(las= 2,
    mar= c(5, 10, 2, 15))
vl_boxplot(diff~res_group, 
           dat, 
           violcol= unique(na.omit(dat[, .(res_group, col)]))$col,
           violin= T,
           ylab= "Obs./Add. (log2)",
           ylim= yl)
abline(h= 0,
       lty= 2)
vl_boxplot(log2FoldChange~res_group, 
           dat, 
           violcol= unique(na.omit(dat[, .(res_group, col)]))$col,
           violin= T,
           ylab= "Activity (log2)")
abline(h= 0,
       lty= 2)
dev.off()
