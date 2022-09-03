setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import motif enrichment
enr <- readRDS("Rdata/vllib002_lm_residuals_SOM_motifs_enr.rds")

pdf("pdf/draft/Figure_2D.pdf", 
    height = 4.5, 
    width = 4)
par(mgp= c(1.5, 0.5, 0),
    mar= c(2,10,1,6),
    tcl= -0.2,
    las= 1)
plot(enr$L,
     padj_cutoff= 0.00001,
     cex.balloons= 0.5, 
     top_enrich= 10, 
     col= c("blue", "red"))
plot(enr$R,
     padj_cutoff= 0.01,
     cex.balloons= 0.5, 
     top_enrich= 6, 
     col= c("blue", "red"))
dev.off()