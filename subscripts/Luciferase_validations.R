setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ----
lib <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- merge(luc,
             lib, 
             by= c("L", "R"),
             suffixes= c("_luc", "_STARR"))
dat[, actClass:= fcase(grepl("^control", L) & grepl("^control", R), "Ctl./Ctl.",
                       grepl("^control", R), "Enh./Ctl.",
                       grepl("^control", L), "Ctl./Enh.",
                       default= "Enh./Enh.")]
dat[, actClass:= factor(actClass, c("Ctl./Ctl.", "Enh./Ctl.", "Ctl./Enh.", "Enh./Enh."))]
dat <- na.omit(dat)
Cc <- c("grey0", "royalblue2", "purple", "#74C27A")
dat[, actCol:= adjustcolor(Cc[actClass], 0.5)]

# Model ----
.lm <- lm(log2FoldChange_luc~log2FoldChange_STARR, dat)

# Plot ----
pdf("pdf/draft/Luciferase_validations.pdf", 
    height = 3, 
    width = 3)
par(mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
dat[, {
  plot(log2FoldChange_STARR, 
       log2FoldChange_luc,
       xaxt= "n",
       xlab= "pSTARR-Seq activity (log2)",
       ylab= "Norm. luciferase activity (log2)",
       ylim= c(-0.8, 6.6),
       col= actCol,
       pch= 16,
       cex= 0.6)
  axis(1, padj = -1.25)
  segments(log2FoldChange_STARR,
           log2FoldChange_luc-sd,
           log2FoldChange_STARR,
           log2FoldChange_luc+sd,
           col= actCol)
  vl_plot_coeff(value = cor.test(log2FoldChange_STARR, 
                                 log2FoldChange_luc)$estimate,
                type = "pcc",
                inset= c(-0.1, 0),
                cex= 7/12)
  .SD
}]
# Legend
abline(.lm, lty= "11")
unique(dat[, .(actClass, actCol)])[,{
  legend("bottomright",
         col= actCol,
         legend= actClass,
         bty= "n",
         cex= 7/12,
         pch= 16,
         y = .75, 
         inset = c(-0.25, 0),
         xpd= NA)
}]
dev.off()