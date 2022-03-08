setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- lib[vllib=="vllib002"]
dat <- feat$add_feature(dat, feat$lib)


pl <- rbind(unique(dat[, .(enh= L, class= class_L, act= median_L, group= group_L, side= "Left")]),
            unique(dat[, .(enh= R, class= class_R, act= median_R, group= group_R, side= "Right")]))
pl[, group:= factor(ifelse(group=="control", "controls", "candidates"), c("controls", "candidates"))]
pl[, class:= fcase(class=="inactive", "inactive", 
                   act<=2, "low [1-2]",
                   act<=4, "medium (2-4]",
                   default= "strong >4")]
pl <- as.matrix(dcast(pl, class~side+group), 1)
Cc <- adjustcolor(c("cornflowerblue", "gold", "tomato", "red"), 0.5)

pdf("pdf/draft/Figure_1D.pdf", width = 4, height = 4.5)
par(las= 2,
    mar= c(7.1, 4.1, 2.1, 6.1), 
    xpd= T)

bar <- barplot(pl,
               col= Cc, 
               ylab= "N oligos",
               ylim= c(0,800),
               names= unlist(tstrsplit(colnames(pl), "_", keep= 2)),
               space= c(1.5,0.1))
seg.y <- par("usr")[1]-strheight("M")*9
segments(bar[c(1, 3)]-strwidth("M")/2, seg.y, bar[c(2, 4)]+strwidth("M")/2, seg.y)
text(c(mean(bar[c(1,2)]), mean(bar[c(3,4)])),
     seg.y,
     pos= 1,
     labels = c("E1", "E2"))
legend(par("usr")[2], 
       par("usr")[4]-strheight("M")*2, 
       fill= rev(Cc),
       legend= rev(rownames(pl)),
       bty= "n",
       xjust= 0,
       cex= 0.8)
text(par("usr")[2],
     par("usr")[4]-strheight("M")*1.5,
     "Individual\nactivity (log2)", 
     pos= 4,
     offset= 0.75)
dev.off()
