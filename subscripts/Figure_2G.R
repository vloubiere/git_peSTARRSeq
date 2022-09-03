setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh.", .(L, R)]
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")[, .(L, R, log2FoldChange, additive)]
dat <- lib[dat, on= c("L", "R"), nomatch= NULL]
dat[, group:= ifelse(grepl("^dev_", L) & grepl("^dev_", R), "Super", "Add")]

pdf("pdf/draft/Figure_2G.pdf", 
    height = 3, 
    width = 3)
par(mgp= c(1.5, 0.5, 0),
    mar= c(3.3,2.5,0.7,7),
    tcl= -0.2,
    las= 1)
# Calls
.b <- substitute(boxplot(x = dat[group=="Super", .(additive, log2FoldChange)], 
                         at= c(1, 2), 
                         boxwex= 0.4,
                         xlim= c(0.5,4), 
                         staplewex= NA, 
                         whisklty= 1, 
                         frame= F,
                         col= adjustcolor("#74C27A", 0.7),
                         pch= 16,
                         outcol= adjustcolor("#74C27A", 0.5),
                         add= T,
                         xaxt= "n"))
.l <- substitute(segments(x0 = 1,
                          y0= dat[group=="Super", additive], 
                          x1 = 2,
                          y1= dat[group=="Super", log2FoldChange],
                          col= adjustcolor("#74C27A", 0.5)))
# Plot
plot.new()
plot.window(xlim= c(0.5, 4.5),
            ylim= c(-0.5, 6.3))
eval(.b)
title(ylab= "Luciferase activity (log2)")
eval(.l)
.b$x <- dat[group=="Add", .(additive, log2FoldChange)]
.b$at <- c(3,4)
.b$col <- adjustcolor("tomato", 0.7)
.b$yaxt <- "n"
eval(.b)
.l$x0 <- 3
.l$y0 <- dat[group=="Add", additive]
.l$x1 <- 4
.l$y1 <- dat[group=="Add", log2FoldChange]
.l$col <- adjustcolor("tomato", 0.5)
eval(.l)
axis(1, at = 1:4, labels = NA)
text(1:4,
     par("usr")[3]-strheight("M", cex= 0.5),
     c("Exp. Add.", "Obs.", "Exp. Add.", "Obs."),
     srt= 45,
     pos= 2,
     xpd= T,
     offset= -0.25)
segments(c(1,3),
         c(6.6, 6.1),
         c(2,4),
         c(6.6, 6.1),
         xpd= T)
pval_sup <- wilcox.test(dat[group=="Super", additive], dat[group=="Super", log2FoldChange])$p.value
pval_add <- wilcox.test(dat[group=="Add", additive], dat[group=="Add", log2FoldChange])$p.value
vl_plot_pval_text(c(1.5, 3.5),
                  c(6.6, 6.1),
                  c(pval_sup, pval_add),
                  stars_only = T,
                  xpd= T)
legend(4,
       par("usr")[4],
       c("Dev. pair", "Hk. contain. pair"),
       fill= adjustcolor(c("#74C27A", "tomato"), 0.7),
       bty= "n",
       xpd= T,
       cex= 0.9)
dev.off()
