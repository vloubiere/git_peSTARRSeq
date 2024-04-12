setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data and match reversed pairs ----
dat <- readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds")
dat <- dat[!grepl("control", L) & !grepl("control", R)]
pair <- merge(dat,
             dat,
             by.x= c("L", "R"),
             by.y= c("R", "L"),
             suffixes= c("", ".r"))

# Collapsed per enhancer
coll <- merge(dat[, .("Mean act. 5' loc. (log2)"= mean(log2FoldChange)), .(enh= L)],
              dat[, .("Mean act. 3' loc. (log2)"= mean(log2FoldChange)), .(enh= R)],
              by= "enh")

# Plot ----
Cc <- c("grey20", "cornflowerblue", "gold", "red")

pdf("pdf/draft/review_orientation_activity.pdf", 6, 3)
vl_par(mfrow= c(1,2),
       pty= "s")
pair[, {
  smoothScatter(log2FoldChange,
                log2FoldChange.r,
                xlab= "A/B activity (log2)",
                ylab= "B/A activity (log2)",
                col= adjustcolor(blues9[9], .3),
                xaxt= "n")
  axis(1, padj= -1.25)
  pcc <- cor.test(log2FoldChange,
                  log2FoldChange.r)$estimate
  vl_plot_coeff(value= pcc,
                type = "pcc",
                cex= 7/12)
  clip(-1.5, 8, -1.5, 8)
  abline(0, 1, lty= "11")
}]
coll[, {
  plot(`Mean act. 3' loc. (log2)`,
       `Mean act. 5' loc. (log2)`,
       col= adjustcolor("lightgrey", .4),
       # col= adjustcolor(Cc[actL], .4),
       pch= 16,
       cex= 0.5,
       xaxt= "n")
  axis(1, padj= -1.25)
  pcc <- cor.test(`Mean act. 3' loc. (log2)`,
                  `Mean act. 5' loc. (log2)`)$estimate
  vl_plot_coeff(value= pcc,
                type = "pcc",
                cex= 7/12)
  clip(1, 7, 1, 7)
  abline(0, 1, lty= "11")
  # legend(par("usr")[2],
  #        par("usr")[4],
  #        cex= 7/12,
  #        bty= "n",
  #        fill= Cc,
  #        legend= levels(actL),
  #        xpd= T,
  #        border= NA)
}]
dev.off()