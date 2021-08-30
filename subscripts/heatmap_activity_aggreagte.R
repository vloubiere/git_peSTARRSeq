setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
dat <- readRDS("Rdata/final_results_table.rds")

# PLOT
pdf("pdf/aggregate_activity_additivity.pdf", height = 4.75, width = 10.5)
par(mfrow= c(1, 2))
dat[, {
  .act <- dcast(data = .SD, 
                formula= group_L~group_R, 
                value.var= "log2FoldChange",
                fun.aggregate= median, na.rm= T)
  .act <- as.matrix(.act, 1)
  .add <- dcast(data = .SD, 
                formula= group_L~group_R, 
                value.var= "diff",
                fun.aggregate= median, na.rm= T)
  .add <- as.matrix(.add, 1)
  row_ord <- order(rowSums(.act))
  col_ord <- order(colSums(.act))
  .act <- .act[row_ord, ]
  .act <- .act[, col_ord]
  .add <- .add[row_ord, ]
  .add <- .add[, col_ord]
  vl_heatmap(mat = .act, 
             cluster_rows = F, 
             cluster_cols = F, 
             col = c("blue", "cornflowerblue", "yellow"),
             main= cdition, 
             display_numbers = T, 
             legend_title = "Activity (log2)")
  vl_heatmap(mat = .add, 
             cluster_rows = F, 
             cluster_cols = F,
             main= cdition, 
             breaks = c(-1.5, -0.5, 0, 0.5, 1.5),
             col = c("cornflowerblue", "white", "white", "white", "tomato"),
             display_numbers = T, 
             legend_title = "Additivity (log2)")
}, cdition]
dev.off()

