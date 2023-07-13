setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
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

pdf("pdf/draft/Luciferase_validations.pdf", 
    height = 3, 
    width = 3)
vl_par(mar= c(3.5, 2.75, 0.5, 1.5),
       mgp= c(1.5, 0.5, 0),
       bty= "n")
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
  .lm <- lm(log2FoldChange_luc~log2FoldChange_STARR, dat)
  abline(.lm, lty= "11")
  PCC <- cor.test(log2FoldChange_STARR, log2FoldChange_luc)$estimate
  leg <- unique(dat[order(actClass), .(actClass, actCol)])
  leg[, legend("bottomright",
               legend = c(paste0("R2= ", round(summary(.lm)$r.squared, 2), " | PCC= ", round(PCC, 2)),
                          as.character(rev(actClass))),
               bty= "n",
               col= c("black", 
                      rev(actCol)),
               pch= c(NA, 
                      rep(16, length(actClass))),
               lty= c("11", 
                      rep(NA, length(actClass))),
               seg.len= 0.5,
               cex= 0.8,
               inset= c(-0.4,0),
               xpd= NA)]
}]
dev.off()

file.show("pdf/draft/Luciferase_validations.pdf")