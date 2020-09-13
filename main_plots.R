setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")

########### Position effect ################
pl <- dat[enh_L != enh_R & !grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R)]
individual <- merge(unique(dat[!is.na(median_L), .(enh= enh_L, median_L)]), unique(dat[!is.na(median_R), .(enh= enh_R, median_R)]), all.x= T, all.y= T)
dir <- lib[individual, median_L:= i.median_L, on= c("ID_vl==enh")]
dir <- dir[individual, median_R:= i.median_R, on= c("ID_vl==enh")]
dir <- dir[!is.na(median_L) & !is.na(median_R)]

pdf("pdf_main/PCC_fwd_rev.pdf", 9, 4)
layout(matrix(1:3, ncol=3), widths = c(1, 1, 0.5))

# All
smoothScatter(pl[, .(log2FoldChange, log2FoldChange_rev)], las= 1, xlab= "X~Y activity (log2)", ylab= "Y~X activity (log2)")
.lm <- lm(log2FoldChange_rev~log2FoldChange, pl)
abline(.lm, lty= 2)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor.test(pl$log2FoldChange, pl$log2FoldChange_rev)$estimate, 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)

# Median L vs median_R
plot(dir[, .(median_L, median_R)], las= 1, xlab= "Left individual activity (log2)", ylab= "Right individual activity (log2)",
     pch= 19, cex= 0.8, col= adjustcolor("grey", 0.6))
.lm <- lm(median_R~median_L, dir)
abline(0, 1)
abline(.lm, lty= 2)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor.test(dir$median_L, dir$median_R)$estimate, 2)
legend("topleft", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n", lty= 2)

# Quantif
dir <- melt(dir, id.vars = "ID_vl", measure.vars = c("median_L", "median_R"))
pirateplot(value~variable, dir, theme= 3, gl.col = NA, ylab= "Activity  (log2)", ylim=c(-2, 10))
my_pval_plot(value~variable, dir)
dev.off()
  
########### additive models ################
# models
sub <- dat[!is.na(log2FoldChange) & !is.na(log2FC_add) & grepl("^dev", enh_L) & grepl("^dev", enh_R)]
linear <- lm(log2FoldChange~log2FC_add, sub)
asymp <- nls(log2FoldChange ~ SSasymp(log2FC_add, Asym, resp0, lrc), data= sub)

# plot
pdf("pdf_main/Asymp_model_obs_vs_exp_add.pdf", 6, 4.5)
layout(matrix(1:2, ncol= 2), widths = c(1, 0.5))
par(mar= c(5.5,4.5,2,1))

smoothScatter(sub[, .(log2FC_add, log2FoldChange)], xlab= "Expected additive activity (log2)", ylab= "Activity (log2)", las=1)
lines(seq(-2, 12, 0.1), predict(asymp, newdata = data.table(log2FC_add= seq(-2, 12, 0.1))), lty= 2)

RMSE_add = sqrt(mean((sub$log2FoldChange-sub$log2FC_add)^2))
RMSE_lm = sqrt(mean((sub$log2FoldChange-predict(linear))^2))
RMSE_asymp = sqrt(mean((sub$log2FoldChange-predict(asymp))^2))

barplot(c(RMSE_add, RMSE_lm, RMSE_asymp), las= 2, ylim= c(0, 1.6), ylab= "Root mean square error", names.arg= c("Additive", "Linear", "Asymptotic"))
dev.off()
  
########### Linear models subgroups ################
sub <- copy(dat)
sub[act_group=="active~active", c("xl", "exp_value") := .("Expected additive activity (log2)", log2FC_add)]
sub[act_group=="active~inactive", c("xl", "exp_value") := .("Median left enhancer", median_L)]
sub[act_group=="inactive~active", c("xl", "exp_value") := .("Median right enhancer", median_R)]
sub <- sub[!is.na(act_group) & !is.na(log2FoldChange) & !is.na(exp_value)]

pdf("pdf_main/Subgroups_linear_models.pdf", 10, 3.35)
layout(matrix(1:4, ncol= 4), widths = c(1,1,1,0.5))
par(mar= c(7,4.5,2,1))
sub[, 
    {
      smoothScatter(exp_value, log2FoldChange, xlab= xl[1], ylab= "Activity (log2)", las= 1, xlim= c(-0, 10), ylim= c(-3, 12), main= act_group[1])
      .lm <- lm(log2FoldChange~exp_value);
      leg <- paste0("PCC= ", round(cor.test(exp_value, log2FoldChange)$estimate, 2), "\n R2= ", round(summary(.lm)$r.squared, 2));
      legend("topleft", legend= leg, bty= "n");
      abline(.lm, lty= 2)
      abline(0, 1)
    }, .(act_group, xl)]
pirateplot(log2FoldChange-exp_value~act_group, sub)

dev.off()
  
########### Heatmap expL expR residuals ################
mat <- as.matrix(dcast(dat[!is.na(diff)], quant_L~quant_R, value.var = "diff", fun.aggregate = median), 1)
mat <- mat[nrow(mat):1,]

pdf("pdf_main/Heatmap_LR_vs_residuals.pdf", 10, height = 9)
layout(matrix(1:4, ncol=2), widths= c(0.25, 1), heights= c(0.25, 1))
par(mar=c(1,1,2.5,1), lty= 2, las= 1)
plot.new()
plot(apply(mat[nrow(mat):1,], 1, mean), seq(nrow(mat)), xlim= rev(c(-0.75, 2)), xaxt= "n", yaxt= "n")
abline(v=0)
axis(3)
par(mar=c(1,2.5,1,6))
plot(apply(mat, 2, mean), xaxt= "n", ylim= c(-0.75, 2))
abline(h=0)
par(mar=c(1.5,3.2,3.2,6.7))
my_pheatmap(mat, cluster_rows= F, cluster_cols=F, lim= c(-2.5, 2.5), col= c("cornflowerblue", "white", "white", "white", "red"),
            legend_title = "log2(o/e)", row_labels = F, col_labels = F, grid_lwd = NA)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], lty= 1)
abline(h=0.0375)
abline(v=0.0375)
dev.off()
  
########### Barplot fraction additive ################
pdf("pdf_main/Density_additive_strong_weak.pdf", 5, 4)
my_empty_plot(xlim= c(-4, 4), ylim= c(0, 0.6), xlab= "observed vs expected additive (log2)", ylab= "density")
strong <- dat[median_L>8 | median_R>8, diff]
lines(density(strong), col= adjustcolor("cornflowerblue", 0.4), lwd= 3)
weak <- dat[median_L>2 & median_L<8 & median_R>2 & median_R<8, diff]
lines(density(weak), col= adjustcolor("tomato", 0.4), lwd= 3)
leg <- c(paste0("Strong pairs (", formatC(length(strong), big.mark = ","), ")"), paste0("Weak pairs (", formatC(length(weak), big.mark = ","), ")"))
legend("topleft", fill= adjustcolor(c("cornflowerblue", "tomato"), 0.4), legend = leg, bty= "n")
dev.off()
  
########### Barplot examples add and sup ################
pl <- dat[ex, .(enh_L, enh_R, median_L, median_R, log2FoldChange), on=c("enh_L", "enh_R")]
pl[feat, tss_L:= i.ctss, on= "enh_L==uniq_ID"]
pl[feat, tss_R:= i.ctss, on= "enh_R==uniq_ID"]
pl[, title:= paste(tss_L, "x", tss_R)]

Cc1 <- c("darkgoldenrod2", "gold1", "midnightblue", "cornflowerblue", "seagreen", "limegreen", "tomato", "lightcoral")
Cc2 <- c("yellow", NA, "mediumturquoise", NA, "greenyellow", NA, "plum1")

pdf("pdf_main/Barplot_examples_peSTARRSeq.pdf", 5.5, 5)
par(mar= c(8,5,5,2))
bar <- barplot(c(2^pl$median_L[1]+2^pl$median_R[1], 2^pl$log2FoldChange[1],
                 2^pl$median_L[2]+2^pl$median_R[2], 2^pl$log2FoldChange[2],
                 2^pl$median_L[3]+2^pl$median_R[3], 2^pl$log2FoldChange[3],
                 2^pl$median_L[4]+2^pl$median_R[4], 2^pl$log2FoldChange[4]),
                 las= 1, ylim= c(0, 800), ylab= "Activity", col = Cc1, space = rep(c(1, 0.1), 4))

text(bar, -30, paste0(rep(pl$title, each= 2), c(" exp.", " obs.")), srt=45, xpd= T, pos= 2, offset= -0.25)

barplot(c(2^pl$median_L[1], NA,
          2^pl$median_L[2], NA,
          2^pl$median_L[3], NA,
          2^pl$median_L[4], NA), add= T, axes= F, col= Cc2, space = rep(c(1, 0.1), 4))

barplot(c(2^pl$median_L[1], NA,
          2^pl$median_L[2], NA,
          2^pl$median_L[3], NA,
          2^pl$median_L[4], NA), add= T, axes= F, density = 10, angle = 45, col= "black", space = rep(c(1, 0.1), 4))
dev.off()
  
########### Scatterplots single enhancers ################
pdf("pdf_main/Scatterplot_sgl_enhancers_LR.pdf", 4, 4.5)
plot(dat[enh_L==SGL, .(median_R, 2^diff)], las= 1, ylab= "observed/expected", xlab= "Activity right candidate (log2)",
     pch= 19, col= adjustcolor("grey", 0.7), main= paste(" sgl enhancer left"), cex= 0.8)
abline(h=1, lty= 2)
plot(dat[enh_R==SGL, .(median_L, 2^diff)], las= 1, ylab= "observed/expected", xlab= "Activity right candidate (log2)",
     pch= 19, col= adjustcolor("grey", 0.7), main= paste("sgl enhancer left"), cex= 0.8)
abline(h=1, lty= 2)
dev.off()
