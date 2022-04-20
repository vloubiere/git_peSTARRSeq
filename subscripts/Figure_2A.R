setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/uniq_enh_feat/lib_genomic_dat.rds")
dat <- lib[vllib=="vllib002" & class== "enh./enh."]
.lm <- lm(log2FoldChange~additive, dat)

pdf("pdf/draft/Figure_2A.pdf", height = 4.5, width = 6)
layout(matrix(1:2, ncol= 2), widths = c(1,0.5))
par(las= 1)
smoothScatter(dat[, .("Expected additive (log2)"= additive,
                      "Activity (log2)"= log2FoldChange)],
              xlim= c(-1.5,9.5),
              ylim= c(-1.5,9.5))
abline(0, 1)
abline(.lm, lty= 2)
text(par("usr")[1], 
     par("usr")[2]-strwidth("M"),
     pos= 4,
     "Super-additive")
text(par("usr")[2], 
     par("usr")[3]+strwidth("M"),
     pos= 2,
     "Sub-additive")
legend(par("usr")[1], 
       par("usr")[2]-strwidth("M")*1.5,
       legend= c("Linear model",
                 paste0("(R²= ", round(summary(.lm)$r.square, 2), ")")),
       lty= c(2, 0),
       bty= "n",
       cex= 0.7)
mtext(paste0("Activity = Expected additive * ", round(.lm$coefficients[2], 1), " + ", round(.lm$coefficients[1], 1)), line=1)
res <- vl_boxplot(list("Exp. add."= dat$additive, 
                       "Observed"= dat$log2FoldChange), 
                  violin= T, 
                  violcol = "lightgrey", 
                  ylab= "Activity (log2)",
                  las= 2,
                  compute_pval = list(c(1,2)))
abline(h= res$stats[3,1], lty= 2)
dev.off()