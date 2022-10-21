setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

dat <- rbindlist(list(DSCP= readRDS("db/FC_tables/vllib015_pe-STARR-Seq-2.0_DSCP_T12_SCR1_300_counts_norm_final_oe.rds"),
                      RpS12= readRDS("db/FC_tables/vllib016_pe-STARR-Seq-2.0_RpS12_T12_SCR1_300_counts_norm_final_oe.rds")), idcol = "CP")
dat <- dat[(grepl("^hk", L) & grepl("^hk", R)) | (grepl("^dev", L) & grepl("^dev", R)) & actClassL=="active" & actClassR=="active"]
dat[, group:= ifelse(grepl("^hk", L), "Hk pairs", "Dev pairs")]

pl <- merge(dat[CP=="DSCP"],
            dat[CP=="RpS12"],
            by= c("L", "R"),
            suffixes= c(".dev", ".hk"))
pl[, Cc:= adjustcolor(ifelse(grepl("hk", L), "tomato", "limegreen"), 0.3)]

pdf("pdf/draft/Compare_add_mult_vllib15_16_activity.pdf",
    height = 3,
    width = 4.25)
par(las= 1,
    mar= c(3,3,1,7),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    cex= 1,
    bty= "n")
pl[, {
  plot(log2FoldChange.dev,
       log2FoldChange.hk,
       col = Cc,
       xlab= "DSCP activity (log2)",
       ylab= "RpS12 activity (log2)",
       pch= 19,
       cex= 0.5)
  abline(0, 1, lty=2)
  PCC <- round(cor.test(log2FoldChange.dev,
                 log2FoldChange.hk)$estimate, 2)
  legend("topleft",
         legend= paste0("PCC= ", PCC),
         bty= "n",
         xpd= T,
         inset= c(-0.1,-0.05))
}]
unique(pl[, .(Cc, group.dev)])[, {
  legend(par("usr")[2],
         par("usr")[4],
         bty= "n",
         pch= 19,
         col= Cc,
         legend= group.dev,
         xpd= T)
}]
par(mfrow= c(1,2),
    mar= c(3,3,1,2))
vl_boxplot(log2FoldChange~CP+group,
           dat,
           names= function(x) gsub(".Dev pairs|.Hk pairs", "", x),
           tilt.names= T,
           compute_pval= list(c(1,2), c(3,4)),
           col= adjustcolor(rep(c("limegreen", "tomato"), each= 2), 0.3),
           ylab= "Activity (log2)",
           at= c(1,2,3.5,4.5))
legend(par("usr")[2],
       par("usr")[4],
       fill= adjustcolor(c("limegreen", "tomato"), 0.3),
       legend= c("Dev/Dev", "Hk/Hk"),
       bty= "n",
       x.intersp = 0.5,
       y.intersp = 0.8,
       cex= 0.8,
       xpd= NA)
dev.off()