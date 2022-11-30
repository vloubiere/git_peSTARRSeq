setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables_DESeq2/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_DESeq2_final_oe.rds")
dat <- dat[actClass %in% c("ctl./ctl.", "enh./ctl.", "ctl./enh.", "enh./enh.")]
dat[, actClass:= droplevels(actClass)]
Cc <- c("grey0", "royalblue2", "purple", "#74C27A")

pdf("pdf/draft/Compare_individual_vs_enh_pairs.pdf", 
    width = 2, 
    height = 3)
par(las= 2, 
    mar= c(3.5,3,0.5,2),
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