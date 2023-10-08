setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import library
dat <- list("Developmental pool (x1,000)"= readRDS("Rdata/vl_library_twist008_112019.rds"),
            "Houskeeping pool (x165)"= readRDS("Rdata/vl_library_twist12_210610.rds")[group %in% c("control", "dev", "hk")])
dat <- lapply(dat, as.data.table)
dat <- rbindlist(dat,
                 idcol = "lib",
                 fill= T)
# Remove Ecoli controls
dat <- dat[!grepl("Ecoli", ID_vl)]
# Compute STARR-Seq signal
dat[, devSTARR:= vl_bw_coverage(dat, "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")]
dat[, hkSTARR:= vl_bw_coverage(dat, "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw")]
# melt and arrange
dat <- melt(dat,
             id.vars = c("lib", "group"),
             measure.vars = patterns("STARR$"))
dat[, group:= factor(group,
                      c("control", "dev", "hk", "ecdysone", "heatshock", "OSC"))]

# Plot
Cc <- adjustcolor(c("limegreen", "tomato"), .5)

pdf("pdf/draft/oligo_pools_characterization.pdf",
    height = 3.5,
    width = 8)
par(las= 1,
    tcl= -.2,
    mgp= c(2, .5, 0),
    mfrow= c(1,2),
    oma= c(0,0,0,10))
dat[, {
  vl_boxplot(log2(value)~variable+group,
             at= rep(seq(1:6), each= 2)+c(-.15, .15),
             xaxt= "n",
             boxwex= 0.15,
             col= Cc,
             ylab= "STARR-Seq signal (log2)",
             main= lib)
  axis(1, labels = NA)
  vl_tilt_xaxis(1:6, 
                labels = levels(group))
  abline(h= 0,
         lty= "11")
  .SD
}, lib]
legend("topright",
       fill= Cc,
       legend= c("Dev. STARR-Seq (DSCP)",
                 "Hk. STARR-Seq (RpS12)"),
       bty= "n",
       cex= .7,
       xpd= NA,
       inset = c(-1, 0))
dev.off()