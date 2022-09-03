setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data
lib <- readRDS("Rdata/final_results_table.rds")[class_act== "enh./enh." & grepl("^hk_|^dev_", L) & grepl("^hk_|^dev_", R)]
lib[, diff:= log2FoldChange-additive]
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
cols <- c("L", "R", "median_L", "median_R", "log2FoldChange", "diff")
dat <- merge(lib[vllib=="vllib015", cols, with= F],
             lib[vllib=="vllib016", cols, with= F],
             by= c("L", "R"),
             suffixes= c("_dev", "_hk"))
dat[grepl("^dev", L) & grepl("^dev", R), c("col", "group"):= .("#74C27A", "Dev./Dev.")]
dat[grepl("^hk", L) & grepl("^hk", R), c("col", "group"):= .("tomato", "Hk./Hk.")]
dat[grepl("^hk", L) & grepl("^dev", R), c("col", "group"):= .("royalblue2", "Hk./Dev.")]
dat[grepl("^dev", L) & grepl("^hk", R), c("col", "group"):= .("goldenrod1", "Dev./Hk.")]
dat <- dat[sample(nrow(dat))]
leg <- unique(dat[, .(group, col)])
leg[, group:= factor(group, c("Dev./Dev.", "Hk./Hk.", "Hk./Dev.", "Dev./Hk."))]
setorderv(leg, "group")
PCC <- dat[, paste0("PCC= ", 
                    round(cor.test(log2FoldChange_dev,
                                   log2FoldChange_hk)$estimate, 2))]
box <- melt(dat, 
            id.vars = "group", 
            measure.vars = c("log2FoldChange_dev", "log2FoldChange_hk"), 
            variable.factor = F)
box <- box[group %in% c("Dev./Dev.", "Hk./Hk.")]
box[, group:= factor(group, c("Dev./Dev.", "Hk./Hk."))]
box[, variable:= switch(variable, 
                        "log2FoldChange_dev"= "dCP",
                        "log2FoldChange_hk"= "hkCP"), variable]
box2 <- melt(dat, 
             id.vars = "group", 
             measure.vars = c("diff_dev", "diff_hk"), 
             variable.factor = F)
box2 <- box2[group %in% c("Dev./Dev.", "Hk./Hk.")]
box2[, group:= factor(group, c("Dev./Dev.", "Hk./Hk."))]
box2[, variable:= switch(variable, 
                         "diff_dev"= "dCP",
                         "diff_hk"= "hkCP"), variable]

pdf("pdf/draft/Figure_3AB.pdf",
    height = 3,
    width= 5)
layout(matrix(1:2, ncol= 2), 
       widths = c(1,.75))
par(mgp= c(1.5, 0.5, 0),
    mar= c(3,3,1,0.5),
    tcl= -0.2,
    las= 1)
dat[, {
  plot(log2FoldChange_dev,
       log2FoldChange_hk,
       col= adjustcolor(col, 0.3),
       pch= 16,
       cex= 0.8,
       xlim= c(-2.5,11.4),
       ylim= c(0,10))
  abline(0,1, lty= 2)
},]
leg[, {
  legend("topleft",
         c(PCC, as.character(group)), 
         col= c(NA, col),
         pch= c(NA, rep(16, .N)),
         bty= "n",
         cex= 0.6)
}]
par(mar= c(3,2.5,1,0.25))
vl_boxplot(value~group+variable, 
           box, 
           tilt.names= T,
           boxcol= rep(c("#74C27A", "tomato"), 2),
           ylab= "Activity (log2)")
legend("topright",
       c("Dev./Dev.", "Hk./Hk."), 
       fill= c("#74C27A", "tomato"),
       bty= "n",
       cex= 0.6)

plot.new()
vl_boxplot(value~group+variable, 
           box2, 
           tilt.names= T,
           boxcol= rep(c("#74C27A", "tomato"), 2),
           ylab= "Observed-exp. additive (log2)")
legend("topright",
       c("Dev./Dev.", "Hk./Hk."), 
       fill= c("#74C27A", "tomato"),
       bty= "n",
       cex= 0.6)
abline(h= 0, lty= 2)
dev.off()

