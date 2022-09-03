setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

cl <- readRDS("Rdata/vllib002_clustering_expected_scores_draft_figure_SOM.rds")

# Check confusion matrix
conf <- merge(cl$rows, cl$cols, by= "name")
tab <- table(conf$cl.x, conf$cl.y)
tab <- matrix(tab, ncol = ncol(tab), dimnames = dimnames(tab))
chi <- chisq.test(tab)
conf <- chi$residuals

# check plot
pdf("pdf/draft/Figure_2C.pdf",
    width = 3.75, 
    height = 3)
par(las= 1, 
    mar= c(2,2,2,6),
    tcl= -0.2,
    mgp= c(2, 0.5, 0))
vl_heatmap(x =  matrix(conf, ncol = ncol(conf), dimnames = dimnames(conf)),
           breaks = c(-2.5,0,2.5),
           cluster_rows = F, 
           cluster_cols = F,
           legend_title = "Standardized\nresidual", 
           display_numbers = T,
           display_numbers_matrix = tab,
           auto_margins = F)
title(main= paste0("Chi-Square pval= ", formatC(chi$p.value, digits = 2)), 
      line= 1,
      font.main= 1,
      cex.main= 0.9)
dev.off()