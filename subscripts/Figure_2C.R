setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
source("git_peSTARRSeq/subscripts/peAddFeatures_function.R")

# Import Clustering object (see clustering_vllib002_lm_residuals_Figure_2B.R)
cl <- readRDS("Rdata/vllib002_clustering_expected_scores_draft_figure_pca.rds")
rows <- cl$rows
cols <- cl$cols
res <- merge(rows, cols, by= "name", suffixes= c("_L", "_R"))
res[, col:= kmeans(as.matrix(res[, .(PC1_L, PC1_R)]), 2)$cluster]
res[, col:= cutree(hclust(dist(as.matrix(SJ(PC1_L, PC1_R)))), 2)]
res[, plot(PC1_L, PC1_R, col= col)]
mot <- fread("Rdata/final_300bp_enhancer_features.txt")
mot <- as.matrix(mot[res$name, grep("motif$", names(mot)), on= "ID", with=F])
enr <- vl_motif_enrich(mot[res$order_L<=50,],
                       mot[res$order_L>=306,], padj_cutoff = 0.01)
enr <- vl_motif_enrich(mot[res$cl_L==1 & res$cl_R==1,],
                       mot[res$cl_L==2 & res$cl_R==2,])

# Contingency table
# tab <- table(paste("Cluster ", res$cl_L), paste0("Cluster ", res$cl_R))
# .f <- fisher.test(tab, alternative = "greater")
# mat <- matrix(c(tab), ncol= 2)
# rownames(mat) <- rownames(tab)
# colnames(mat) <- colnames(tab)
# vl_heatmap(mat, 
#            cluster_rows= F, 
#            cluster_cols= F, 
#            show_col_clusters = F,
#            show_row_clusters = F, 
#            display_numbers = T, 
#            display_numbers_cex = 1, 
#            col = c("cornflowerblue", "tomato"), 
#            legend_title = "N candidates")
# box()
# abline(v=1.5)
# abline(h=1.5)
# title(main= paste0("OR= ", round(.f$estimate, 1), " (pval= ", formatC(.f$p.value, digits = 1), "; Fisher.test)"))

# PCA pcc
pdf("test/test.pdf")
par(mfrow=c (2,2))
# res[order(motif), {
res[, {
  # Cc <- colorRampPalette(c("lightgrey", "gold", "tomato", "red"))(max(motif)+1)
  plot(PC1_L, 
       PC1_R, 
       pch= 16, 
       cex= 1.5,
       # col = adjustcolor(Cc[motif+1], 0.6))
       # col = Cc[motif+1])
       col = kmeans_col)
  legend("topleft",
         paste0("PCC= ", round(cor.test(PC1_L, PC1_R)$estimate, 2)),
         bty= "n")
  abline(h= 0, lty= "11")
  abline(v= 0, lty= "11")
  abline(0, 1, lty= "11")
}]

dev.off()
