setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds")
dat[, actClass:= fcase(grepl("^control", L) & grepl("^control", R), "Ctl./Ctl.",
                       grepl("^control", R), "Enh./Ctl.",
                       grepl("^control", L), "Ctl./Enh.",
                       default= "Enh./Enh.")]
dat[, actClass:= factor(actClass, c("Ctl./Ctl.", "Enh./Ctl.", "Ctl./Enh.", "Enh./Enh."))]
Cc <- c("grey0", "royalblue2", "purple", "#74C27A")

# Plot ----
pdf("pdf/draft/Compare_individual_vs_enh_pairs.pdf", 
    width = 3, 
    height = 3)
# Scatter plot
par(mai= c(0.75,1.25,1,1.1), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75)
vl_boxplot(log2FoldChange~actClass,
           dat,
           compute.pval = list(c(1,2), c(1,3), c(2, 3), c(2,4), c(3,4)),
           pval.cex= 7/12,
           col = adjustcolor(Cc, 0.5),
           ylab= "Activity (log2)",
           tilt.names= T,
           notch= T)
abline(h= 0, 
       lty= "11")
dev.off()
