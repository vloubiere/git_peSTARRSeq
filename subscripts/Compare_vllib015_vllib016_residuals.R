setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

dat <- rbindlist(list(DSCP= readRDS("db/FC_tables/vllib015_pe-STARR-Seq-2.0_DSCP_T12_SCR1_300_counts_norm_final_oe.rds"),
                      RpS12= readRDS("db/FC_tables/vllib016_pe-STARR-Seq-2.0_RpS12_T12_SCR1_300_counts_norm_final_oe.rds")), idcol = "CP")
dat <- dat[(grepl("^hk", L) & grepl("^hk", R)) | (grepl("^dev", L) & grepl("^dev", R)) & actClassL=="active" & actClassR=="active"]
dat[, group:= ifelse(grepl("^hk", L), "Hk", "Dev")]

pl <- melt(dat, 
           id.vars = c("CP", "group", "log2FoldChange"), 
           measure.vars = c("additive", "multiplicative"))
pl[, residuals:= log2FoldChange-value]

pdf("pdf/draft/Compare_add_mult_vllib15_16.pdf",
    height = 3.25,
    width = 4)
par(las= 1,
    oma= c(0,0,0,4),
    mar= c(4.5,3,1.5,0.5),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    mfrow= c(1,2),
    cex= 1)
pl[, {
  vl_boxplot(residuals~variable+group, 
             .SD,
             names= function(x) gsub(".Dev|.Hk", "", x),
             tilt.names= T,
             compute_pval= list(c(1,2), c(3,4)),
             col= adjustcolor(rep(c("limegreen", "tomato"), each= 2), 0.3),
             ylab= "Observerd/Expected (log2)",
             at= c(1,2,3.5,4.5),
             main= CP)
  abline(h= 0,
         lty= 2)
  if(.GRP==.NGRP)
    legend(par("usr")[2],
           par("usr")[4],
           fill= adjustcolor(c("limegreen", "tomato"), 0.3),
           legend= c("Dev/Dev", "Hk/Hk"),
           bty= "n",
           x.intersp = 0.5,
           y.intersp = 0.8,
           cex= 0.8,
           xpd= NA)
  .SD
}, CP]
dev.off()
