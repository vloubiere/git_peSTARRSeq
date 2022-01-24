setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
dat <- readRDS("Rdata/final_results_table.rds")
dat <- dat[L!=R]
# Filter for minimum N
dat[, check:= .SD[, rep(.N>100, .N), .(group_L, group_R)]$V1, vllib]
dat <- dat[(check), !"check"]

# PLOT
pdf("pdf/analyses/balloons_plots_aggregate_activity_additivity.pdf",
    height = 4.5,
    width = 14.5)
par(mfrow= c(1, 3))
dat[, {
  # Activity
  .act <- dcast(data = .SD, 
                formula= group_L~group_R, 
                value.var= "log2FoldChange",
                fun.aggregate= median, na.rm= T)
  .diff <- dcast(data = .SD, 
                 formula= group_L~group_R, 
                 value.var= "diff",
                 fun.aggregate= median, na.rm= T)
  vl_balloons_plot(as.matrix(.act, 1),
                   as.matrix(.diff, 1),
                   x_breaks= c(-2, 7),
                   color_breaks= c(-2, 0, 2),
                   main= paste(c(as.character(vllib), CP, spacer), collapse = " "), 
                   balloon_size_legend= "Activity",
                   balloon_col_legend= "o/e")
}, keyby= .(vllib, CP, spacer)]
dev.off()
