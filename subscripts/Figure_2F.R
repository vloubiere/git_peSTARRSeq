setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(randomForest)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
cl <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")
lib <- readRDS("Rdata/final_results_table.rds")
feat <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
enr <- readRDS("Rdata/Enriched_motifs_high_low_synergy_Figure_2.rds")

# Merge luc and STARR-Seq 
dat <- lib[vllib=="vllib002" & class=="enh./enh.", .(L,R)][luc, on= c("L", "R"), nomatch= NULL]
dat[cl$rows, `5' cluster`:= cl, on= "L==name"]
dat[cl$cols, `3' cluster`:= cl, on= "R==name"]
dat[, col:= fcase(`5' cluster`==1 & `3' cluster`==1, "cyan",
                  `5' cluster`==1 & `3' cluster`==2, "cornflowerblue",
                  # `5' cluster`==2 & `3' cluster`==1, "royalblue2", Does not exist in the data!
                  `5' cluster`==2 & `3' cluster`==2, "tomato")]
dat[, group:= ifelse(col=="tomato", "5'cl.2/3'cl.2", "5'cl.1/3'cl.1\n5'cl.1/3'cl.2")]
dat[, group:= factor(group, c("5'cl.2/3'cl.2", "5'cl.1/3'cl.1\n5'cl.1/3'cl.2"))]
dat[, col:= adjustcolor(col, 0.7)]
dat <- dat[order(!`5' cluster`==2 & `3' cluster`==2)]

pdf("pdf/draft/Figure_2F.pdf", 
    height = 4.5, 
    width = 3)
x <- list("Exp. add."= dat[group=="5'cl.2/3'cl.2", additive],
          "Observed"= dat[group=="5'cl.2/3'cl.2", log2FoldChange],
          "Exp. add."= dat[group!="5'cl.2/3'cl.2", additive],
          "Observed"= dat[group!="5'cl.2/3'cl.2", log2FoldChange])
vl_boxplot(x,
           compute_pval = list(c(1,2), c(3,4)),
           ylab= "Luciferase activity (log2)",
           boxwex = 0.3, 
           las= 1,
           ylim= c(3,6.5),
           xlim= c(0.75, 4.25))
set.seed(2)
jit <- jitter(rep(0, nrow(dat)), factor = 4)
segments(jit+rep(c(1,3), lengths(x[c(1,3)])), 
         unlist(x[c(1,3)]), 
         jit+rep(c(2,4), lengths(x[c(2,4)])),
         unlist(x[c(2,4)]),
         col= "grey40", 
         lwd= 0.25)
points(jit+rep(c(1,3), lengths(x[c(1,3)])), 
       unlist(x[c(1,3)]),
       col= adjustcolor("lightgrey", 0.7), 
       pch= 16,
       cex= 0.5)
points(jit+rep(c(2,4), lengths(x[c(2,4)])), 
       unlist(x[c(2,4)]),
       col= adjustcolor("lightgrey", 0.7), 
       pch= 16,
       cex= 0.6)
dev.off()