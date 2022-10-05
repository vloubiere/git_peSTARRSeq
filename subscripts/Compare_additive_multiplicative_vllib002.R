setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & !grepl("^control", L) & !grepl("^control", R)]

pdf("pdf/draft/Compare_add_mult_vllib002.pdf",
    height = 3,
    width = 1.5)
par(las= 1,
    mar= c(4,3,0,0.5),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    cex= 1)
lib[, {
  vl_boxplot(log2FoldChange-additive,
             log2FoldChange-multiplicative,
             tilt.names= T,
             ylab= "Observed/Expected (log2)",
             names= c("Additive", "Multiplicative"), 
             notch= T,
             col= "grey", 
             compute_pval= list(c(1,2)))
  abline(h= 0, lty= "11")
}]
dev.off()
