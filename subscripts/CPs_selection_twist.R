dat <- as.data.table(readRDS("/groups/stark/lorbeer/vincent/STAP_oligo.RDS"))

#------------------#
# Dev
#------------------#
# plot(log2(dat[["no_enh_Rep1"]]+1), log2(dat[["zfh1_enh_avg"]]+1))
# sel_dev <- identify(log2(dat[["no_enh_Rep1"]]+1), log2(dat[["zfh1_enh_avg"]]+1))
sel_dev <- c(2, 981, 3162, 3857, 3947, 4819, 8818, 9375, 11774, 26652, 18923, 24334, 
             25708, 31105, 41280, 44252, 48942, 51347, 54557, 58040, 61217, 64029, 
             67145, 69447, 11860, 18002, 36510, 51113, 57730, 68507)
plot(log2(dat[["no_enh_Rep1"]]+1), log2(dat[["zfh1_enh_avg"]]+1), 
     col= ifelse(seq(nrow(dat)) %in% sel_dev, "red", adjustcolor("lightgrey", 0.3)), 
     pch= 19)
text(log2(dat[sel_dev][["no_enh_Rep1"]]+1), log2(dat[sel_dev][["zfh1_enh_avg"]]+1))

#------------------#
# Hk
#------------------#
# plot(log2(dat[["no_enh_Rep1"]]+1), log2(dat[["ssp3_enh_Rep1"]]+1))
# sel_hk <- identify(log2(dat[["no_enh_Rep1"]]+1), log2(dat[["ssp3_enh_Rep1"]]+1))
sel_hk <- c(7073, 11537, 15507, 17537, 23126, 24299, 24875, 47428, 48975, 54656, 61815, 62736, 63998, 68989, 69747, 70358)
plot(log2(dat[["no_enh_Rep1"]]+1), log2(dat[["ssp3_enh_Rep1"]]+1), 
     col= ifelse(seq(nrow(dat)) %in% sel_hk, "red", adjustcolor("lightgrey", 0.3)), 
     pch= 19)
text(log2(dat[sel_hk][["no_enh_Rep1"]]+1), log2(dat[sel_hk][["ssp3_enh_Rep1"]]+1))

saveRDS(dat[c(sel_dev, sel_hk)], "Rdata/CPs_selection_twist.rds")
