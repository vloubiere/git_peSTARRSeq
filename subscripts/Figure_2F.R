setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import
pca <- readRDS("Rdata/vllib002_lm_residuals_PCA.rds")
res <- merge(pca$rows,
             pca$cols,
             by= "name", 
             suffixes= c("_L", "_R"))
feat <- fread("Rdata/final_300bp_enhancer_features.txt")[, c("ID", "group", "col")]
res[feat, col:= i.col, on= "name==ID"]
leg <- unique(feat[ID %in% res$name, .(group, col)])

pdf("pdf/draft/Figure_2F.pdf",
    width= 3,
    height= 3)
par(mgp= c(1.5, 0.5, 0),
    mar= c(3,3,1,0.75),
    tcl= -0.2,
    las= 1)
res[, {
  plot(PC1_R, 
       PC1_L,
       pch= 16,
       xlab= "3' PC1",
       ylab= "5' PC1",
       col= adjustcolor(col, 0.6))
  abline(0, 1, lty = 3)
  legend("bottomright",
         c(paste0("PCC= ", round(cor.test(PC1_R, PC1_L)$estimate, 2)), 
           leg$group),
         lty= c(3, rep(0, nrow(leg))),
         pch= c(NA, rep(16, nrow(leg))),
         bty= "n",
         cex= 0.6,
         col= c("black", adjustcolor(leg$col, 0.6)),
         pt.cex= 1,
         seg.len= 1)
}]
dev.off()