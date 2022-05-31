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
                   act<=2, "low (<=2)",
                   act<=4, "medium (2-4)",
                   default= "strong (>4)")]
PCC <- dcast(pl, enh+group~side, value.var = "act")
PCC[dat, TWIST:= dev_log2FC_TWIST_L, on= "enh==L"]
PCC[dat, TWIST_col:= fcase(TWIST<=2 | is.na(TWIST), "#4D9221",
                           TWIST<=6, "#B8E186",
                           TWIST<=8, "#F1B6DA",
                           TWIST>8, "#C51B7D")]
PCC[group=="controls", TWIST_col:= "grey20"]
setorderv(PCC, "group", -1)
.lm <- lm(PCC[, .(Left, Right)])
categ <- as.matrix(dcast(pl, class~side+group, fun.aggregate = length), 1)
Cc <- adjustcolor(c("cornflowerblue", "gold", "tomato", "red"), 0.5)

#------------------------------------------------------#
# PLOT
#------------------------------------------------------#
pdf("pdf/draft/Figure_1CD.pdf", 
    width = 5.75, 
    height = 3)
# Scatter plot
par(mar= c(3.5, 3, 0.5, 0.25), 
    mfrow= c(1,2),
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2)
plot(PCC[, .(Right, Left)],
     pch= 16,
     cex= 0.5,
     col= adjustcolor(PCC$TWIST_col, 0.7),
     las= 1,
     xlab= "3' individual activity (log2)",
     ylab= "5' individual activity (log2)",
     xaxt= "n",
     yaxt= "n")
axis(1, 
     lwd= 0,
     lwd.ticks= 1)
axis(2, 
     lwd= 0,
     lwd.ticks= 1)
text(par("usr")[1],
     par("usr")[4]-strheight("M")/1.5,
     "Candidates",
     pos= 4, 
     cex= 0.7,
     offset= 0.25)
leg <- legend("topleft",
              legend = c("",
                         "strong (>8)",
                         "medium (6-8)",
                         "weak (2-6)",
                         "inact. (<=2)/NA"),
              pch= c(NA, 16, 16, 16, 16, 16),
              col= adjustcolor(c(NA,"#C51B7D","#F1B6DA","#B8E186","#4D9221"), 0.8),
              box.lty= "11",
              cex= 0.5,
              seg.len= 1)
legend(leg$rect$left+leg$rect$w,
       leg$rect$top,
       legend = "control",
       pch= 16,
       col= adjustcolor("grey20", 0.8),
       bty= "n",
       cex= 0.5,
       seg.len= 1)
legend("bottomright",
       legend = paste0("R²= ", round(summary(.lm)$r.squared, 2)," (PCC= ", round(cor.test(PCC$Left, PCC$Right)$estimate, 2), ")"),
       lty= 1,
       bty= "n",
       cex= 0.6,
       seg.len= 1)
abline(h= 0, lty= "11")
segments(0,
         par("usr")[1], 
         0,
         leg$rect$top-leg$rect$h, 
         lty= "11")
abline(.lm)

# Barplot
par(mar= c(3.5, 4, 0.5, 5),
    mgp= c(2, 0.5, 0))
bar <- barplot(categ,
               col= Cc, 
               ylab= "N oligos",
               ylim= c(0,800),
               space= c(1,0.1),
               las= 2,
               border= NA,
               xaxt= "n")
text(bar,
     par("usr")[3]-strheight("M")/2,
     pos= 2,
     unlist(tstrsplit(colnames(categ), "_", keep= 2)),
     xpd= T,
     srt= 45,
     offset= -0.15)
seg.y <- par("usr")[1]-strheight("M")*9
segments(grconvertX(0.15556, "nfc", "user"), 
         grconvertY(0.06833, "nfc", "user"),
         grconvertX(0.2155, "nfc", "user"), 
         grconvertY(0.01333, "nfc", "user"),
         xpd= T)
segments(grconvertX(0.15556, "nfc", "user")+diff(bar[c(1,3)]), 
         grconvertY(0.06833, "nfc", "user"),
         grconvertX(0.2155, "nfc", "user")+diff(bar[c(1,3)]), 
         grconvertY(0.01333, "nfc", "user"),
         xpd= T)
text(grconvertX(0.15556, "nfc", "user")+c(0, diff(bar[c(1,3)])),
     mean(grconvertY(0.01333, "nfc", "user")),
     labels = c("5'", "3'"),
     xpd= T)
leg <- legend(par("usr")[2], 
              par("usr")[4]-strheight("M", cex= 0.8)*2, 
              fill= rev(Cc),
              legend= rev(rownames(categ)),
              bty= "n",
              xjust= 0,
              cex= 0.7,
              xpd= T,
              border= NA)
text(leg$rect$left,
     par("usr")[4]-strheight("M", cex= 0.9),
     "Individual\nactivity (log2)", 
     pos= 4,
     offset= 0.5,
     xpd= T,
     cex= 0.8)
dev.off()
# file.show("pdf/draft/Figure_1CD.pdf")