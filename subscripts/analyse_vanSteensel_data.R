dat <- fread("db/public/GSE240586_EEP_normalized_data_MMA2023030806.tsv.gz")


# Compute my scores
setnames(dat,
         c("id1", "id2"),
         c("L", "R"))
# pseudo <- min(dat$activity_all[dat$activity_all>0])
dat[, log2Foldchange:= log2(activity_all+0.2)]
dat[, indP:= mean(log2Foldchange[suffix1=="Controls" & suffix2=="Controls"], na.rm= T), Promoter]
dat[, indL:= mean(log2Foldchange[suffix2=="Controls"]), .(L, Promoter)]
dat[, indR:= mean(log2Foldchange[suffix1=="Controls"]), .(R, Promoter)]

pdf("test/VanSteensel.pdf", 3, 11)
par(mfcol= c(8, 2),
    mar= c(3,3,1,1),
    bty= "n",
    las= 1,
    tcl= -0.2,
    cex= 0.8,
    cex.axis= 0.8,
    cex.lab= 0.8,
    cex.main= 0.6,
    mgp= c(1.5, 0.15, 0))
dat[, {
  smoothScatter(dat[, log2(2^indL+2^indR-2^indP)],
                dat$log2Foldchange,
                ylab= "Additive",
                xlab= "Observed",
                main= Promoter)
  abline(0, 1, lty= "11")
}, Promoter]
dat[, {
  smoothScatter(dat[, log2(2^indL*2^indR/2^indP)],
                dat$log2Foldchange,
                ylab= "Multiplicative",
                xlab= "Observed",
                main= Promoter)
  abline(0, 1, lty= "11")
}, Promoter]
dev.off()