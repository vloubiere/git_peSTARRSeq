setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
cl <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")
lib <- readRDS("Rdata/final_results_table.rds")
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")

# Merge luc and STARR-Seq 
dat <- lib[vllib=="vllib002" & class=="enh./enh.", .(L,R)][luc, on= c("L", "R"), nomatch= NULL]
dat[cl$rows, `5' cluster`:= cl, on= "L==name"]
dat[cl$cols, `3' cluster`:= cl, on= "R==name"]

pdf("pdf/draft/Figure_2D.pdf", 
    height = 3, 
    width = 1.5)
par(mgp= c(1.5, 0.5, 0),
    mar= c(3.3,2.5,0.7,0.25),
    tcl= -0.2,
    las= 1)
x <- list("Exp. add."= dat[`5' cluster`=="A" & `5' cluster`=="A", additive],
          "Observed"= dat[`5' cluster`=="A" & `5' cluster`=="A", log2FoldChange],
          "Exp. add."= dat[!(`5' cluster`=="A" & `5' cluster`=="A"), additive],
          "Observed"= dat[!(`5' cluster`=="A" & `5' cluster`=="A"), log2FoldChange])
box <- vl_boxplot(x,
                  compute_pval = list(c(1,2), c(3,4)),
                  ylab= "Luciferase activity (log2)",
                  boxwex = 0.3, 
                  las= 1,
                  xlim= c(0.75, 4.25), 
                  tilt.names = T)
set.seed(2)
jit <- jitter(rep(0, nrow(dat)), factor = 4)
segments(rep(c(1.25,3.25), lengths(x[c(1,3)])), 
         unlist(x[c(1,3)]), 
         rep(c(1.75,3.75), lengths(x[c(2,4)])),
         unlist(x[c(2,4)]),
         col= "grey40", 
         lwd= 0.5)
points(rep(c(1.25,3.25), lengths(x[c(1,3)])), 
       unlist(x[c(1,3)]),
       col= adjustcolor("lightgrey", 0.9), 
       pch= 16,
       cex= 0.3)
points(rep(c(1.75,3.75), lengths(x[c(2,4)])), 
       unlist(x[c(2,4)]),
       col= adjustcolor("lightgrey", 0.9), 
       pch= 16,
       cex= 0.3)
text(c(1.5, 3.5), 
     box$pval$y, 
     c("A/A", "Others"),
     pos= 3,
     xpd= T)
dev.off()