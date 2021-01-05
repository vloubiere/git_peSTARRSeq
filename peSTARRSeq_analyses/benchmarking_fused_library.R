setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

#----------------------------------------------------------------#
# 1- Import data 
#----------------------------------------------------------------#
lib <- as.data.table(readRDS("Rdata/library/uniq_library_final.rds"))
lib <- lib[(vl)]
if(!exists(".c"))
{
  all_counts <- readRDS("db/read_counts/all_uniq_counts.rds")
  all_counts[, simp_sample:= ifelse(grepl("^input", sample), "input", "DSCP")]

  .c <- all_counts[grepl("input", sample)]
  .c[lib, group_L:= i.group, on= "enh_L==ID"]
  .c[lib, group_R:= i.group, on= "enh_R==ID"]
  .c[lib, sub_L:= i.linker_ID, on= "enh_L==ID"]
  .c[lib, sub_R:= i.linker_ID, on= "enh_R==ID"]
}

saturation <- CJ(enh_L= lib$ID, enh_R= lib$ID)
.m <- all_counts[simp_sample=="input", sum(counts), .(enh_L, enh_R)]
saturation[.m, input:= i.V1, on= c("enh_L", "enh_R")]
.m <- all_counts[simp_sample=="DSCP", sum(counts), .(enh_L, enh_R)]
saturation[.m, DSCP:= i.V1, on= c("enh_L", "enh_R")]
saturation[is.na(input), input:= 0]
saturation[is.na(DSCP), DSCP:= 0]

dat <- readRDS("Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2.rds")
dat[lib, group_L:= i.group, on= "enh_L==ID"]
dat[lib, group_R:= i.group, on= "enh_R==ID"]

#----------------------------------------------------------------#
# 2- Fusion bias?
#----------------------------------------------------------------#
# sub-library
rdm_sub <- CJ(enh_L= lib$ID, enh_R= lib$ID)
rdm_sub[lib, sub_L:= i.linker_ID, on= "enh_L==ID"]
rdm_sub[lib, sub_R:= i.linker_ID, on= "enh_R==ID"]
rdm_sub <- rdm_sub[, .(random_frac= .N/nrow(rdm_sub)), .(sub_L, sub_R)]

sub <- na.omit(.c)
sub <- sub[, .(frac= sum(counts)/sum(sub$counts)), .(sub_L, sub_R)]
sub[rdm_sub, random_frac:= i.random_frac, on= c("sub_L", "sub_R")]
sub[, ratio := log2(frac)-log2(random_frac)]
sub <- dcast(sub, sub_L~sub_R, value.var = "ratio")

# Activity/enhancer type
rdm_act <- CJ(enh_L= lib$ID, enh_R= lib$ID)
rdm_act[lib, group_L:= i.group, on= "enh_L==ID"]
rdm_act[lib, group_R:= i.group, on= "enh_R==ID"]
rdm_act <- rdm_act[, .(random_frac= .N/nrow(rdm_act)), .(group_L, group_R)]

act <- na.omit(.c)
act <- act[, .(frac= sum(counts)/sum(act$counts)), .(group_L, group_R)]
act[rdm_act, random_frac:= i.random_frac, on= c("group_L", "group_R")]
act[, ratio := log2(frac)-log2(random_frac)]
dact <- dcast(act, group_L~group_R, value.var = "ratio")

# Composition
.N <- na.omit(dat)
.cN <- .N[, .N, .(group_L, group_R)]
dN <- dcast(.cN, group_L~group_R, value.var = "N")

# Simplified version of groups
.sN <- copy(.N)
.sN[, group_L:= switch(group_L, "control"= "control", "dev"= "developmental", "hk"= "housekeeping", "inducible"= "inducible", "OSC"= "non-cell-type-specific", "shared"= "developmental"), group_L]
.sN[, group_R:= switch(group_R, "control"= "control", "dev"= "developmental", "hk"= "housekeeping", "inducible"= "inducible", "OSC"= "non-cell-type-specific", "shared"= "developmental"), group_R]
.sN <- .sN[, .N, .(group_L, group_R)]
dsN <- dcast(.sN, group_L~group_R, value.var = "N")

#----------------------------------------------------------------#
# 3- PLOT
#----------------------------------------------------------------#
pdf("pdf/peSTARRSeq/fused_library_biases.pdf", 9, 12)
par(mar= c(5,5,3,7), mfrow= c(3,2))

# sub libraries
my_pheatmap(as.matrix(sub, 1), cluster_rows = F, cluster_cols = F, display_numbers = T, 
            col= c("cornflowerblue", "white", "tomato"), breaks = c(-1,0,1), legend_title = "log2(o/e)", legend_cex = 0.8)
mtext("sub-library bias", line = 1)
my_fig_label("A", cex= 2)

# Activity
my_pheatmap(as.matrix(dact, 1), cluster_rows = F, cluster_cols = F, display_numbers = T, 
            col= c("cornflowerblue", "white", "tomato"), breaks = c(-1 , 0, 1), legend_title = "log2(o/e)", legend_cex = 0.8)
mtext("candidate type/activity bias", line = 1)
my_fig_label("B", cex= 2)

# Saturation
plot(NA, pch= NA, xlim = c(0, 1e6), ylim = c(0,15), las=1, xlab= "combinations", ylab= "log2(counts+1)")
lines(sort(log2(saturation$input+1)), lwd=3, col= adjustcolor("cornflowerblue", 0.6))
lines(sort(log2(saturation$DSCP+1)), lwd=3, col= adjustcolor("tomato", 0.6))
us <- paste("N=", formatC(nrow(dat), big.mark= ",", format = "fg"), "usable pairs")
legend("topleft", lwd= c(3,3,NA), col = adjustcolor(c("cornflowerblue", "tomato"), 0.6), 
       legend= c("input", "DSCP", us), bty= "n")
my_fig_label("C", cex= 2)

# Number of pairs
my_pheatmap(as.matrix(dN, 1), cluster_rows = F, cluster_cols = F, display_numbers = T, 
            col= c("cornflowerblue", "yellow"), legend_title = "Number", ROUND_FUN = function(x){formatC(x, big.mark = ",")},
            legend_cex = 0.8)
mtext("candidate type combinations", line = 1)
my_fig_label("D", cex= 2)

# Number of pairs
par(mar= c(10,10,3,7))
my_pheatmap(as.matrix(dsN, 1), cluster_rows = F, cluster_cols = F, display_numbers = T, 
            col= c("cornflowerblue", "yellow"), legend_title = "Number", ROUND_FUN = function(x){formatC(x, big.mark = ",")},
            legend_cex = 0.8)
mtext("candidate type combinations", line = 1)
my_fig_label("D", cex= 2)

# Saturation
par(mar= c(5,5,3,7))
plot(NA, pch= NA, xlim = c(-0.5, 12), ylim = c(0,0.6), las=1, xlab= "log2(counts+1)", ylab= "density")
lines(density(log2(saturation$input+1)), lwd=3, col= adjustcolor("cornflowerblue", 0.6))
lines(density(log2(saturation$DSCP+1)), lwd=3, col= adjustcolor("tomato", 0.6))
legend("topright", lwd= c(3,3), col = adjustcolor(c("cornflowerblue", "tomato"), 0.6), 
       legend= c("input", "DSCP"), bty= "n")
my_fig_label("E", cex= 2)

dev.off()
