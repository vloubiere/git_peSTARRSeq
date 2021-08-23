setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#-------------------------------#
# Compute FC
#-------------------------------#

dat <- data.table(file= list.files("db/merged_counts", "vllib016", full.names = T))
dat <- dat[, fread(file), file]
dat[, cdition:= tstrsplit(basename(file), "peSTARRSeq_|_merged", keep= 2), file]
dat <- dat[type=="pair"]

mat <- dcast(dat,
             L+R~cdition,
             value.var = "umi_counts",
             fill = 0)
mat <- mat[, .(L, R, DSCP= DSCP_rep1+DSCP_rep2, input= input_rep1)]
mat <- mat[rowSums(mat[, DSCP:input])>10]
cols <- grep("DSCP|input", names(mat), value = T)
mat[, (cols):= lapply(.SD, function(x) (x+1)/sum(x)*1e6), .SDcols= cols]
mat[, FC:= DSCP/input]

res <- merge(mat[!xor(grepl("^control", L), grepl("^control", R)), .(L, R, FC)],
             mat[grepl("^control", R), .(median_L= median(FC)), L],
             by= "L")
res <- merge(res,
             mat[grepl("^control", L), .(median_R= median(FC)), R],
             by= "R")

#-------------------------------#
# plots
#-------------------------------#
# ALL
smoothScatter(log2(res$median_L+res$median_R),
              log2(res$FC))
abline(0, 1)
sub <- res[grepl("dev", L) & grepl("dev", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("hk", L) & grepl("hk", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("CP", L) & grepl("CP", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[(grepl("dev", L) & grepl("CP", R)) | (grepl("CP", L) & grepl("dev", R))]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("CP", L) & grepl("CP", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[(grepl("dev", L) & grepl("Repressor", R)) | (grepl("Repressor", L) & grepl("dev", R))]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[(grepl("dev", L) & grepl("DHS_peak", R)) | (grepl("DHS_peak", L) & grepl("dev", R))]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[(grepl("dev", L) & grepl("hk", R)) | (grepl("hk", L) & grepl("dev", R))]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("hk", L) & grepl("hk", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("DHS_peak", L) & grepl("DHS_peak", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("DHS_peak", L) & grepl("CP", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("dev", L) & grepl("SUHW", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("dev", R) & grepl("SUHW", L)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
sub <- res[grepl("Repressor", L) & grepl("Repressor", R)]
smoothScatter(log2(sub$median_L+sub$median_R),
              log2(sub$FC))
abline(0, 1)
