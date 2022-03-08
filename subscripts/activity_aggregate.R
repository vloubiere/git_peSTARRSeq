setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
dat <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat <- dat[L!=R]
dat <- feat$add_feature(dat, feat$lib)
dat[, diff:= log2FoldChange-log2(2^(median_L)+2^(median_R))]
# Filter for minimum N
dat[, check:= .SD[, rep(.N>100, .N), .(group_L, group_R)]$V1, vllib]
dat <- dat[(check), !"check"]

# PLOT
pdf("pdf/analyses/balloons_plots_aggregate_activity_additivity.pdf",
    height = 4.5,
    width = 15.5)
par(mfrow= c(1, 3),
    las= 2)
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
                   color_breaks= c(-2, 0, 2),
                   x_breaks= c(-2, 0, 2, 4, 6),
                   main= paste(c(as.character(vllib), CP, spacer), collapse = " "), 
                   balloon_size_legend= "Activity",
                   balloon_col_legend= "o/e", 
                   cex.balloons= 0.8)
}, keyby= .(vllib, CP, spacer)]
dev.off()
