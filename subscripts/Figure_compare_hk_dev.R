setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(vioplot)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
lib[, diff:= log2FoldChange-additive]
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")

# dat
dat <- merge(lib[vllib=="vllib015"][act_wilcox_L<0.05 & act_wilcox_R<0.05, 
                                    .(L, R, log2FoldChange, diff)],
             lib[vllib=="vllib016"][act_wilcox_L<0.05 & act_wilcox_R<0.05, 
                                    .(L, R, log2FoldChange, diff)],
             suffixes= c("_dev", "_hk"),
             by= c("L", "R"))
# Compute colors and legend
group_cols <- feat$add_feature(dat, feat$lib)
group_cols[, col:= fcase(group_L=="dev" & group_R=="dev", col_L,
                         group_L=="dev" & group_R=="hk", "gold",
                         group_L=="hk" & group_R=="dev", "royalblue2",
                         group_L=="hk" & group_R=="hk", col_L,
                         default= "lightgrey")]
dat[group_cols, col:= col,on= c("L", "R")]
leg <- unique(group_cols[order(col=="lightgrey"), .(paste0(group_L," x ",group_R), col)][, .(name= ifelse(col=="lightgrey", "Others", V1), col)])
# Format before plotting
dat <- melt(dat,
            measure.vars = c("log2FoldChange_hk", "log2FoldChange_dev", "diff_hk", "diff_dev"))
dat[, c("variable", "group"):= tstrsplit(variable, "_")]
dat <- dcast(dat, L+R+col+variable~group)
dat <- dat[order(col!="lightgrey")]

# plot
pdf("pdf/analyses/DSCP_RpS12_compare_act_residuals.pdf", 5, 5)
par(mar= c(5,5,3,3),
    las= 1, 
    lty= 2,
    bty= "n")
dmat[, {
  # Scatter
  plot(dev, 
       hk, 
       col= adjustcolor(col, 0.3),
       xlab= switch(variable,
                    "log2FoldChange"= "Dev. Activity (log2)",
                    "diff"= "Dev. Observed/Additive (log2)"),
       ylab= switch(variable,
                    "log2FoldChange"= "Hk. Activity (log2)",
                    "diff"= "Hk. Observed/Additive (log2)"),
       pch= 16)
  abline(0,1)
  abline(h= 0)
  abline(v= 0)
  
  # Densities
  bottom <- par("usr")[4]
  top <- grconvertY(1, "nfc", "user")-strheight("M")/2
  dev_dens <- .SD[col!="lightgrey", density(dev)[c("x", "y")], col]
  dev_dens <- dev_dens[order(factor(col, levels= c("gold", "royalblue2", "tomato", "#74C27A")))]
  dev_dens[, adj.y:= bottom+y/max(y)*(top-bottom), col]
  dev_dens[, polygon(x, adj.y, xpd= T, border=NA, col= adjustcolor(col[1], 0.5)), col]
  
  left <- par("usr")[2]
  right <- grconvertX(1, "nfc", "user")-strwidth("M")/2
  hk_dens <- .SD[col!="lightgrey", setNames(density(hk)[c("x", "y")], c("y", "x")), col]
  hk_dens <- hk_dens[order(factor(col, levels= c("gold", "royalblue2", "tomato", "#74C27A")))]
  hk_dens[, adj.x:= left+x/max(x)*(right-left), col]
  hk_dens[, polygon(adj.x, y, xpd= T, border=NA, col= adjustcolor(col[1], 0.5)), col]
  
  # Legends
  legend("topleft", 
         bty= "n",
         fill= leg$col,
         legend = leg$name,
         border= NA,
         cex= 0.75)
  par(bty= "o")
  box()
}, variable]
dev.off()

