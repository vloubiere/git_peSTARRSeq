setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- lib[vllib=="vllib002"]
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- merge(luc,
             dat, 
             by= c("L", "R"),
             suffixes= c("_luc", "_STARR"),
             all.x= T)
dat <- na.omit(dat)
dat[, "class":= fcase(class %in% c(NA, "ctl./ctl.", "inact./ctl.", "ctl./inact.", "inact./inact."), "inact./inact.",
                      class %in% c("ctl./enh.", "inact./enh."), "inact./enh.",
                      class %in% c("enh./ctl.", "enh./inact."), "enh./inact.",
                      class=="enh./enh.", "enh./enh.")]
dat[, class:= factor(class, c("inact./inact.", "enh./inact.", "inact./enh.", "enh./enh."))]
Cc <- c(adjustcolor("lightgrey", 0.5), "#436EEE80", adjustcolor("purple", 0.5), "#74C27A80")
dat[, col:= Cc[.GRP], keyby= class]
dat <- dat[!is.na(log2FoldChange_luc)]
dat[, class:= droplevels(class)]

leg <- unique(dat[order(class), .(class, col)])
.lm <- lm(log2FoldChange_STARR~log2FoldChange_luc, dat)

pdf("pdf/draft/Figure_1FG.pdf", 
    height = 4.5, 
    width = 7.2)
layout(matrix(1:2, ncol= 2), widths = c(1, 0.65))
par(las= 1,
    mar= c(5.1, 4.1, 2.1, 2.1))
dat[, {
  plot(log2FoldChange_STARR, 
       log2FoldChange_luc,
       xlab= "pe-STARR-Seq activity (log2)",
       ylab= "Normalized luciferase activity (log2)",
       ylim= c(-0.8, 6.6),
       col= col,
       pch= 16)
  segments(log2FoldChange_STARR,
           log2FoldChange_luc-sd,
           log2FoldChange_STARR,
           log2FoldChange_luc+sd,
           col= col)
}]
abline(.lm, lty=2)
legend("topleft",
       legend = c(paste0("R²= ", round(summary(.lm)$r.squared, 2), 
                         " (PCC= ", round(cor.test(dat$log2FoldChange_STARR, dat$log2FoldChange_luc)$estimate, 2), ")"),
                  levels(leg$class)),
       bty= "n",
       col= c("black", Cc),
       lty= c(2,0,0,0,0),
       pch= c(NA,19,19,19,19),
       cex= 0.65)
vl_boxplot(log2FoldChange_luc~class,
           dat, 
           violin= T,
           las= 2, 
           compute_pval= list(c(1,2), c(2,4), c(3,4)),
           violcol= Cc,
           ylab= "Normalized luciferase activity (log2)", 
           ylab.line = 2)
dev.off()
