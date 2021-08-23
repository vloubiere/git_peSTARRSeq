setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#------------------------#
# Import 
#------------------------#
lib <- readRDS("Rdata/final_300bp_enhancer_features.rds")
all <- fread("db/final_tables_exp_model/vllib013_short_spacer_peSTARRSeq_max_cutoff_final_oe.txt")
dat <- na.omit(all[median_L>0.5 & median_R>0.5 & L != R & !grepl("^control", L) & !grepl("^control", R)])
cols <- c("ID", grep("^motif__", colnames(lib), value = T))
dat <- rbind(merge(dat, 
                   lib[, ..cols], 
                   by.x= "L",
                   by.y= "ID"),
             merge(dat, 
                   lib[, ..cols], 
                   by.x= "R",
                   by.y= "ID"))
dat <- dat[, lapply(.SD, function(x) log2(sum(x)+1)), L:median_R]

#------------------------#
# Build models 
#------------------------#
#------- Additive -------#
dat[, exp:= log2(2^median_L+2^median_R)]
.lm1 <- lm(log2FoldChange~exp, 
           dat)
dat[, pred1:= predict(.lm1)]
#------- Linear model -------#
.lm2 <- lm(log2FoldChange~median_L*median_R, 
           dat)
dat[, pred2:= predict(.lm2)]
#------- With motifs -------#
.flm3 <- paste0(c("log2FoldChange~median_L*median_R", grep("^motif", colnames(dat), value= T)), collapse= "+")
.lm3 <- lm(as.formula(.flm3), 
           dat)
dat[, pred3:= predict(.lm3)]
fwrite(dat, "Rdata/vllib013_modelling_sum_motif_counts.txt")

#------------------------#
# lm3 coefficients for barplot
#------------------------#
best_fit <- as.data.table(summary(.lm3)$coefficients, keep.rownames = T)
best_fit[, FDR:= p.adjust(`Pr(>|t|)`, method = "fdr")]
best_fit <- best_fit[FDR<0.05][order(`t value`)]
best_fit[, rn:= gsub("median_L", "A", rn)]
best_fit[, rn:= gsub("median_R", "B", rn)]
best_fit[, name:= gsub("^motif__", "", rn)]
som <- readRDS("Rdata/som_clustering_motifs_300bp_enhancers.rds")$info
setkeyv(som, "motif")
best_fit[name %in% som$motif, name:= som[name, BA_cluster], name]

#------------------------#
# Melt for boxplot
#------------------------#
pl <- melt(dat[, .(exp, log2FoldChange)], 
           measure.vars =  c("exp", "log2FoldChange"))
pl[, variable:= as.character(variable)]
pl[, variable:= switch(variable, exp= "Add. (log2)", log2FoldChange= "Obs. (log2)"), variable]

#------------------------------------------------------------------------#
# PLOT
#------------------------------------------------------------------------#
dir.create("pdf/modeling/", showWarnings = F)

pdf("pdf/modeling/vllib013_modeling_sum_motif_counts.pdf", 
    width = 15.5, 
    height = 8.5)

lims <- c(-1, 11)
layout(matrix(c(1:5, rep(6, 5)), ncol= 5, nrow = 2, byrow = T), 
       widths = c(1, 0.4, 1, 1, 1))

#------ MODEL 1 ------#
par(las= 1,
    mar= c(8,6,4,0.5))
smoothScatter(dat[, .(exp, log2FoldChange)], 
              xlab= "Additive (A+B, log2)",
              ylab= "Observed (AB, log2)",
              xlim= lims,
              ylim= lims)
mtext("(A+B)", 
      line = 1.5)
.f1 <- gsub("log2FoldChange", "AB", vl_model_equation(.lm1, 
                                                      digits= 2))
.f1 <- gsub("exp", "(A+B)", .f1)
mtext(.f1, 
      line = 0.1,
      cex= 0.7)
abline(0, 1)
abline(.lm1, lty= 2)
legend("topleft", 
       legend = c(paste0("PCC= ", round(cor.test(dat$log2FoldChange, dat$exp)$estimate, 2)), 
                  paste0("R2= ", round(summary(.lm1)$r.squared, 2))),
       bty= "n")

#------ BOXPLOT ------#
par(las= 2,
    mar= c(8,4,4,0.5))
vl_boxplot(x = na.omit(pl), 
           formula= value~variable, 
           col= "black", 
           ylab = "Activity (log2)")

#------ MODEL 2 ------#
par(las= 1,
    mar= c(8,6,4,0.5))
smoothScatter(dat[, .(pred2, log2FoldChange)], 
              xlab= "Predicted (log2)",
              ylab= "Observed (log2)",
              xlim= lims,
              ylim= lims)
mtext("A*B", 
      line = 1.5)
.f2 <- gsub("log2FoldChange", "AB", vl_model_equation(.lm2, 
                                                      digits= 2))
.f2 <- gsub("median_L", "A", .f2)
.f2 <- gsub("median_R", "B", .f2)
mtext(.f2, 
      line = 0.1,
      cex= 0.7)
abline(0, 1, lty= 2)
legend("topleft", 
       legend = c(paste0("PCC= ", round(cor.test(dat$log2FoldChange, dat$pred2)$estimate, 2)), 
                  paste0("R2= ", round(summary(.lm2)$r.squared, 2))),
       bty= "n")

#------ MODEL 3 ------#
smoothScatter(dat[, .(pred3, log2FoldChange)], 
              xlab= "Predicted (log2)",
              ylab= "Observed (log2)",
              xlim= lims,
              ylim= lims)
mtext("A*B+mot1+mot2+...", 
      line = 1.5)
.f3 <- gsub("log2FoldChange", "AB", vl_model_equation(.lm3, 
                                                      digits= 2))
.f3 <- gsub("motif__", "", .f3)
.f3 <- substr(.f3, 1, 80)
.f3 <- gsub(".flm3", "AB", .f3)
.f3 <- gsub("median_L", "A", .f3)
.f3 <- gsub("median_R", "B", .f3)
mtext(paste0(.f3, "..."), 
      line = 0.1, 
      cex= 0.7)
abline(0, 1, lty= 2)
legend("topleft", 
       legend = c(paste0("PCC= ", round(cor.test(dat$log2FoldChange, dat$pred3)$estimate, 2)), 
                  paste0("R2= ", round(summary(.lm3)$r.squared, 2))),
       bty= "n")

#------ barplot coefficients ------#
par(las= 2)
barplot(best_fit[, `t value`], 
        names.arg = best_fit$name, 
        ylim= c(-20, 60),
        ylab= "Linear model t value", 
        border= NA,
        col = ifelse(best_fit[, `t value`]>0, "tomato", "cornflowerblue"))
dev.off()

