setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(kohonen)


cl <- readRDS("Rdata/clustering_lm_residuals_vllib002.rds")

###############################################
# PLOT 
###############################################
confusion <- table(merge(cl$rows, cl$cols, by= "name")[, .(cl.x, cl.y)])
chi <- chisq.test(confusion)
chi.res <- chi$residuals
mat <- matrix(chi.res, nrow(chi.res), ncol(chi.res), dimnames = dimnames(chi.res))

pdf("pdf/draft/confusion_matrix_L_R_clusters.pdf", 4, 3)
par(mar= c(6,6.75,1,5.25),
    tcl= -0.2,
    las= 2)
vl_heatmap(mat, 
           cluster_rows = F, 
           cluster_cols = F, 
           display_numbers = T,
           display_numbers_matrix = matrix(confusion, nrow(confusion), ncol(confusion)),
           breaks= c(-4.5, 0, 4.5), 
           legend_title = "Standardized\nresidual")
dev.off()