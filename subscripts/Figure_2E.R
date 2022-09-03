setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import Clustering object and retrieve enhacner sequences
cl <- readRDS("Rdata/vllib002_lm_residuals_SOM.rds")
rows <- cl$rows
cols <- cl$cols
feat <- fread("Rdata/final_300bp_enhancer_features.txt")[, .(ID, ATAC, K27Ac, K4me1, K4me3)]
rows <- feat[rows, on= "ID==name"]
cols <- feat[cols, on= "ID==name"]

pdf("pdf/draft/Figure_2E.pdf",
    width= 7,
    height= 3)
par(mfrow=c(1,4),
    mgp= c(3, 0.5, 0),
    mar= c(2,5,2,1),
    tcl= -0.2,
    las= 1)
rows[, {
  lapply(seq(.SD), function(i)
  {
    vl_boxplot(.SD[[i]]~as.character(cl),
               rows, 
               compute_pval= list(c(1,2), c(1,3), c(1,4)),
               ylab= "Enrichment",
               boxcol= c("lightgrey", "darkgrey", "tomato", "cornflowerblue")[i],
               main= c("ATAC-Seq", "H3K27Ac", "H3K4me1", "H3K4me3")[i])
  })
}, .SDcols= c("ATAC", "K27Ac", "K4me1", "K4me3")]
dev.off()