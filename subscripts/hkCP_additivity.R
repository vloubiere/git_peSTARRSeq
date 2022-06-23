setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features_wo_motifs.rds")
dat <- lib[vllib=="vllib016" & class== "enh./enh."]
dat[feat, group_L:= i.group, on= "L==ID"]
dat[feat, group_R:= i.group, on= "R==ID"]
dat[group_L=="dev" & group_R=="dev", c("col", "class"):= .("#74C27A", "dev/dev")]
dat[group_L=="dev" & group_R=="hk", c("col", "class"):= .("cyan", "dev/hk")]
dat[group_L=="hk" & group_R=="dev", c("col", "class"):= .("royalblue2", "hk/dev")]
dat[group_L=="hk" & group_R=="hk", c("col", "class"):= .("tomato", "hk/hk")]
dat[, class:= factor(class, c("dev/dev", "hk/dev", "dev/hk", "hk/hk"))]
dat <- dat[!is.na(class)]

pdf("pdf/draft/hkCP_activity.pdf", 
    height = 3, 
    width = 3.8)
layout(matrix(1:2, ncol= 2), 
       widths = c(1,0.35))
par(las= 1,
    mar= c(3.5,2.75,0.5,0.1),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2)
plot(NA, 
     xlim= c(-1, 10.2),
     ylim= c(-1, 10.2),
     xlab= "Expected Additive (log2)",
     ylab= "Activity (log2)")
dat[, points(additive, 
             log2FoldChange,
             pch= 16,
             col= adjustcolor(col, 0.5),
             cex= 0.25)]
dat[, {
  biv_kde <- MASS::kde2d(additive,
                         log2FoldChange,
                         h= 2)
  contour(biv_kde,
          col= colorspace::darken(col[1], amount = 0.3),
          add = T,
          drawlabels = F,
          nlevels = 5,
          lwd= 0.5)
}, .(class, col)]

# smoothScatter(dat[, .("Expected additive (log2)"= additive,
#                       "Activity (log2)"= log2FoldChange)],
#               xlim= c(-1.5,12),
#               ylim= c(-1.5,12),
#               colramp = colorRampPalette(c("white", "grey90", "grey60", "grey20")),
#               xaxt= "n",
#               yaxt= "n")
# axis(1, 
#      lwd= 0,
#      lwd.ticks= 1)
# axis(2, 
#      lwd=0,
#      lwd.ticks=1)
abline(0, 1)
# abline(.lm, lty= 2)
# text(par("usr")[1], 
#      par("usr")[2]-strwidth("M"),
#      pos= 4,
#      "Super-additive")
# text(par("usr")[2], 
#      par("usr")[3]+strwidth("M"),
#      pos= 2,
#      "Sub-additive")
# legend(par("usr")[1],
#        par("usr")[2]-strwidth("M")*0.9,
#        legend= c("Linear model",
#                  paste0("(R²= ", round(summary(.lm)$r.square, 2), ")")),
#        lty= c(2, 0),
#        bty= "n",
#        cex= 0.7,
#        x.intersp= 0.5)
# mtext(paste0("Activity = Expected additive * ", round(.lm$coefficients[2], 1), " + ", round(.lm$coefficients[1], 1)), line=1)
# 
# # boxplot
# par(xpd= T)
# res <- vl_boxplot(list("Exp. add."= dat$additive, 
#                        "Observed"= dat$log2FoldChange), 
#                   violin= T, 
#                   violcol = "lightgrey", 
#                   ylab= "Activity (log2)",
#                   las= 2,
#                   compute_pval = list(c(1,2)),
#                   ylim= c(-3, 13),
#                   tilt.names = T)
# abline(h= res$stats[3,1], 
#        lty= 2,
#        xpd= F)
dev.off()






test <- dcast(dat[group_L=="hk"], median_L+L~median_R+R, value.var = "log2FoldChange")
vl_heatmap(as.matrix(test[, .(median_L, `0.70005383848753_dev_A_01303`, `2.60554939811674_hk_A_01333`, `6.47809379076528_hk_A_01362`)]),
           cluster_rows=F, 
           cluster_cols= F)




