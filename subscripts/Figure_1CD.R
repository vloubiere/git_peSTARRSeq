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
                   act<=2, "low",
                   act<=4, "medium",
                   default= "strong")]
PCC <- dcast(pl, enh+group~side, value.var = "act")
setorderv(PCC, "group", -1)
.lm <- lm(PCC[, .(Left, Right)])
categ <- as.matrix(dcast(pl, class~side+group), 1)
Cc <- adjustcolor(c("cornflowerblue", "gold", "tomato", "red"), 0.5)

#------------------------------------------------------#
# PLOT
#------------------------------------------------------#
pdf("pdf/draft/Figure_1CD.pdf", width = 8, height = 5)
# Scatter plot
par(mar= c(7.1, 4.1, 2.1, 0.1), 
    mfrow= c(1,2))
plot(PCC[, .(Right, Left)],
     pch= 16,
     col= adjustcolor(PCC[, ifelse(group=="candidates", "grey80", "grey20")], 0.5),
     las= 1,
     xlab= "3' individual activity (log2)",
     ylab= "5' individual activity (log2)")
abline(v= 1, lty= 2)
abline(h= 1, lty= 2)
legend("topleft",
       legend = c("control", "candidate"),
       pch= 16,
       col= c("grey20", "grey80"),
       bty= "n")
abline(.lm)
legend("bottomright",
       legend = paste0("R²= ", round(summary(.lm)$r.squared, 2)," (PCC= ", round(cor.test(PCC$Left, PCC$Right)$estimate, 2), ")"),
       lty= c(1,0),
       bty= "n",
       cex= 0.8)

# Barplot
par(mar= c(7.1, 6.1, 2.1, 6.1))
bar <- barplot(categ,
               col= Cc, 
               ylab= "N oligos",
               ylim= c(0,800),
               names= unlist(tstrsplit(colnames(categ), "_", keep= 2)),
               space= c(1.5,0.1),
               las= 2)
seg.y <- par("usr")[1]-strheight("M")*9
segments(bar[c(1, 3)]-strwidth("M")/2, 
         seg.y, 
         bar[c(2, 4)]+strwidth("M")/2, 
         seg.y,
         xpd= T)
text(c(mean(bar[c(1,2)]), mean(bar[c(3,4)])),
     seg.y,
     pos= 1,
     labels = c("5'", "3'"),
     xpd= T)
leg <- legend(par("usr")[2], 
              par("usr")[4]-strheight("M")*2, 
              fill= rev(Cc),
              legend= rev(rownames(categ)),
              bty= "n",
              xjust= 0,
              cex= 0.8,
              xpd= T)
text(leg$rect$left,
     par("usr")[4]-strheight("M")*1.5,
     "Individual\nactivity (log2)", 
     pos= 4,
     offset= 0.5,
     xpd= T)
dev.off()
