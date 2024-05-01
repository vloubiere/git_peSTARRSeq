setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
S2 <- readRDS("db/FC_tables/DSCP_large_WT_FC_DESeq2.rds")
ECD <- readRDS("db/FC_tables/DSCP_ECD_WT_FC_DESeq2.rds")
OSC <- readRDS("db/FC_tables/DSCP_OSC_WT_FC_DESeq2.rds")

# Select pairs ----
pl <- list("ECD - ECD"= S2[grepl("^ecd", L) & grepl("^ecd", R), log2FoldChange],
           "ECD + ECD"= ECD[grepl("^ecd", L) & grepl("^ecd", R), log2FoldChange],
           "OSC in S2"= S2[grepl("^OSC", L) & grepl("^OSC", R), log2FoldChange],
           "OSC in OSC"= OSC[grepl("^OSC", L) & grepl("^OSC", R), log2FoldChange])

# Melt ----
Cc <- c("lightgrey", "cornflowerblue", "lightgrey", "limegreen")
pdf("pdf/draft/ecd_OSC_induction.pdf", 3, 3)
vl_par(mgp= c(1, .35, 0))
vl_boxplot(pl,
           at= c(0.9,1.1,1.9,2.1),
           xaxt= "n",
           ylab= "Activity (log2)",
           lwd= .75,
           boxwex= .2,
           col= adjustcolor(Cc, .6))
abline(h= 0,
       lty= "13",
       lwd= .75)
vl_tilt_xaxis(x = c(1,2),
              labels= c(paste0("Ecd. inducible enh.\npairs (n= ", length(pl$`ECD + ECD`), ")"),
                        paste0("OSC-specific enh.\npairs (n= ", length(pl$`OSC in OSC`), ")")))
legend(par("usr")[2]-strwidth("M"),
       par("usr")[4],
       fill= adjustcolor(Cc[-3], .6),
       legend = c("S2 cells - ecdysone",
                  "S2 cells + ecdysone",
                  "OSC cells"),
       xpd= T,
       bty= "n",
       cex= 7/12,
       border= NA)
dev.off()