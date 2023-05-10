setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
dat[, actClass:= fcase(grepl("^control", L) & grepl("^control", R), "ctl./ctl.",
                       actL!="Inactive" & actR!="Inactive", "enh./enh.",
                       grepl("^control", L) & actR!="Inactive", "ctl./enh.",
                       actL!="Inactive" & grepl("^control", R), "enh./ctl.", 
                       default = NA)]
dat[, actClass:= factor(actClass, c("ctl./ctl.", "enh./ctl.", "ctl./enh.", "enh./enh."))]
Cc <- c("grey0", "royalblue2", "purple", "#74C27A")

pdf("pdf/draft/Compare_individual_vs_enh_pairs.pdf", 
    width = 2, 
    height = 3)
par(las= 2, 
    mar= c(3.5,3,0.75,2),
    tcl= -0.2,
    mgp= c(1.5, 0.5, 0),
    lty= 1)
dat[, {
  vl_boxplot(log2FoldChange~actClass,
             compute_pval = list(c(1,2), c(1,3), c(2,4), c(3,4)),
             col = adjustcolor(Cc, 0.5),
             ylab= "Activity (log2)",
             tilt.names= T,
             notch= T)
  print("")
}]
abline(h= 0, 
       lty= 2)
dev.off()

file.show("pdf/draft/Compare_individual_vs_enh_pairs.pdf")