require(data.table)
require(GenomicRanges)
require(rtracklayer)
require(yarrr)
setwd("~/Dropbox (VBC)/R_db/dm3/0002_library_design/")

dCP <- import("../Available_data/STARR-Seq/peaks/dCP_S2_merged_peaks.bed")
dist <- distanceToNearest(dCP)
dist <- data.table("dCP", as.data.table(dist)$distance)

pirateplot(V2~V1, dist, ylim = c(0, 100000))
