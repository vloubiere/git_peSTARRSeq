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

# cuts <- seq(-1, 9, 0.1)
# grid <- CJ(x0= cuts[-length(cuts)],
#            y0= cuts[-length(cuts)])
# grid[, x1:= x0+0.1]
# grid[, y1:= y0+0.1]
# grid[, x:= rowMeans(.SD), .SDcols= c("x0", "x1")]
# grid[, y:= rowMeans(.SD), .SDcols= c("y0", "y1")]
# grid$N <- dat[grid, .N, .EACHI, on= c("additive>x0", "additive<=x1", "log2FoldChange>y0", "log2FoldChange<=y1")]$N
# mat <- dcast(grid, y~x, value.var = "N")
# mat <- as.matrix(mat, 1)
# Cc <- circlize::colorRamp2(c(0, max(mat)),
#                            c("white", "black"))
# im <- Cc(mat)
# par(las= 1)
# plot.new()
# plot.window(xlim= c(-1, 9),
#             ylim= c(-1, 9))
# rasterImage(im,
#             xleft = 0, 
#             ybottom = 9,
#             xright= 9, 
#             ytop= 0)

# diags <- col(mat)+row(mat)
# diag_idx <- unique(c(diags))
# bg_cols <- colorRampPalette(c("cornflowerblue", "white", "tomato"))
# for(i in diags)
#   diags[diags==i] <- adjustcolor(bg_cols(max(diag_idx))[i], 0.5)
# diags <- diags[, ncol(diags):1]
# rasterImage(diags,
#             xleft = 0, 
#             ybottom = 9,
#             xright= 9, 
#             ytop= 0)

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