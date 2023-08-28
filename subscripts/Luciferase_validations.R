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
dat[, actClass:= fcase(grepl("^control", L) & grepl("^control", R), "ctl./ctl.",
                       actL!="Inactive" & actR!="Inactive", "enh./enh.",
                       grepl("^control", L) & actR!="Inactive", "ctl./enh.",
                       actL!="Inactive" & grepl("^control", R), "enh./ctl.",
                       default = NA)]
dat[, actClass:= factor(actClass, c("ctl./ctl.", "enh./ctl.", "ctl./enh.", "enh./enh."))]
dat <- na.omit(dat)
Cc <- c("grey0", "royalblue2", "purple", "#74C27A")
dat[, actCol:= adjustcolor(Cc[actClass], 0.5)]

# Model ----
.lm <- lm(log2FoldChange_luc~log2FoldChange_STARR, dat)

# Plot ----
pdf("pdf/draft/Luciferase_validations.pdf", 
    height = 3, 
    width = 3)
par(mar= c(3.5, 2.75, 0.5, 1.5),
    mgp= c(1.5, 0.5, 0),
    bty= "n",
    las= 1,
    tcl= -0.2)
dat[, {
  plot(log2FoldChange_STARR, 
       log2FoldChange_luc,
       xlab= "pSTARR-Seq activity (log2)",
       ylab= "Norm. luciferase activity (log2)",
       ylim= c(-0.8, 6.6),
       col= actCol,
       pch= 16,
       cex= 0.8)
  segments(log2FoldChange_STARR,
           log2FoldChange_luc-sd,
           log2FoldChange_STARR,
           log2FoldChange_luc+sd,
           col= actCol)
  .SD
}]
# Legend
abline(.lm, lty= "11")
unique(dat[, .(actClass, actCol)])[,{
  legend("bottomright",
         fill= actCol,
         legend= actClass,
         bty= "n",
         cex= 0.8)
}]
vl_plot_R2(rsquare = summary(.lm)$r.squared,
           inset= c(-0.1,0))
dev.off()