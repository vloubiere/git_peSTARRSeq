setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table.rds")
dat[, diff:= log2(2^median_L+2^median_R)]
dat <- dat[grepl("^hk|^dev", L) & grepl("^hk|^dev", R)]
unique(dat[, .(vllib, library,spacer, CP)])

dev <- merge(dat[vllib=="vllib023", !c("vllib","library","CP","spacer","class","col")], # short intron 4
             dat[vllib=="vllib018", !c("vllib","library","CP","spacer","class","col")], # 2kb intron 4
             c("L", "R"),
             suffixes= c("_s", "_l"))
hk <- merge(dat[vllib=="vllib024", !c("vllib","library","CP","spacer","class","col")], # short intron 4
            dat[vllib=="vllib020", !c("vllib","library","CP","spacer","class","col")], # 2kb intron 4
            c("L", "R"),
            suffixes= c("", "_i")) 

par(mfrow= c(2,1))
smoothScatter(dev$log2FoldChange_s, dev$log2FoldChange_l)
abline(0, 1, lty= 2)
smoothScatter(dev$diff_s, dev$diff_l)
abline(0, 1, lty= 2)

smoothScatter(hk$log2FoldChange, hk$log2FoldChange_i)
abline(0, 1, lty= 2)


mdev <- merge(unique(dev[, .(enh= L, median_L, median_L_i, act_wilcox_L_i<0.05)]),
              unique(dev[, .(enh= R, median_R, median_R_i, act_wilcox_R_i<0.05)]),
              by= "enh")
mhk <- merge(unique(hk[, .(enh= L, median_L, median_L_i, act_wilcox_L_i<0.05)]),
             unique(hk[, .(enh= R, median_R, median_R_i, act_wilcox_R_i<0.05)]),
             by= "enh")
vl_boxplot(list(mdev[,median_R-median_L],
                mdev[,median_R_i-median_L_i]), 
           compute_pval = list(c(1,2)),
           notch= T)
abline(h= 0)
vl_boxplot(list(mhk[,median_R-median_L],
                mhk[,median_R_i-median_L_i]), 
           compute_pval = list(c(1,2)),
           notch= T)
abline(h= 0)

plot(NA, 
      xlim=c(-4,10),
      ylim=c(-4,10))
points(dat[vllib=="vllib018" & class== "enh./enh.", .(additive, log2FoldChange)], col= adjustcolor("green", 0.5))
abline(0, 1)

plot(NA, 
     xlim=c(-4,10),
     ylim=c(-4,10))
points(dat[vllib=="vllib020" & class== "enh./enh.", .(additive, log2FoldChange)], col= adjustcolor(dat[vllib=="vllib020" & class== "enh./enh.", col], 0.5))
abline(0, 1)
