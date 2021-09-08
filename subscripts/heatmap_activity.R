setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")

dat <- DT[L!=R]
dat <- dat[!is.na(active_plot_group_LR)]
dat <- dat[!is.na(additive) & !is.na(log2FoldChange) & !is.na(diff)]

# PLOT
pdf("pdf/all_pairs_activity_additivity.pdf", height = 4.75, width = 5*3)
par(mfrow= c(1, 3))
dat[, {
  # Use expected additive values to order the heatmap
  .add <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "additive")
  .add <- as.matrix(.add, 1)
  row_ord <- order(-rowMeans(.add, na.rm = T))
  col_ord <- order(colMeans(.add, na.rm = T))
  .add <- .add[row_ord,]
  .add <- .add[, col_ord]
  vl_heatmap(mat = .add, 
             main= active_plot_group_LR, 
             cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             show_colnames = F,
             col = viridis::plasma(10),
             legend_title = "Exp. Add. (log2)")
  rect(0,0,1,1,lwd= 0.5)
  # Activity
  .act <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "log2FoldChange")
  .act <- as.matrix(.act, 1)
  .act <- .act[row_ord,]
  .act <- .act[, col_ord]
  vl_heatmap(mat = .act, 
             main= active_plot_group_LR, 
             cluster_rows = F,
             cluster_cols = F,
             col = viridis::plasma(10),
             show_rownames = F,
             show_colnames = F,
             legend_title = "Activity (log2)")
  rect(0,0,1,1,lwd= 0.5)
  # Additivity
  .diff <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "diff")
  .diff <- as.matrix(.diff, 1)
  .diff <- .diff[row_ord,]
  .diff <- .diff[, col_ord]
  .breaks <- max(abs(quantile(.diff, c(0.05, 0.975), na.rm= T)))
  vl_heatmap(mat = .diff, 
             main= active_plot_group_LR, 
             cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             show_colnames = F,
             breaks = c(-.breaks, -0.5, 0.5, .breaks),
             col= c("cornflowerblue", "white", "white", "tomato"),
             legend_title = "Additivity (log2)")
  rect(0,0,1,1,lwd= 0.5)
  # # Motif
  # .mot <- dcast(data = .SD[, .(L, R, value= log2(motif__flyfactorsurvey__suHw_FlyReg_FBgn0003567_R+1))], 
  #                formula= L~R, 
  #                value.var= "value")
  # .mot <- as.matrix(.mot, 1)
  # .mot <- .mot[row_ord,]
  # .mot <- .mot[, col_ord]
  # .breaks <- max(abs(quantile(.mot, c(0.05, 0.975), na.rm= T)))
  # vl_heatmap(mat = .mot, 
  #            main= plot_group, 
  #            cluster_rows = F,
  #            cluster_cols = F,
  #            show_rownames = F,
  #            show_colnames = F,
  #            # col= c("cornflowerblue", "white", "white", "tomato"),
  #            col = viridis::plasma(10),
  #            # breaks = c(-.breaks, -0.5, 0.5, .breaks),
  #            legend_title = "Motif counts")
  # rect(0,0,1,1)
  print("DONE")
}, active_plot_group_LR]
dev.off()


