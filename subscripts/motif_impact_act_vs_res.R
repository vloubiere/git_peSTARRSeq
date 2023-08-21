setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
dat$ctlL <- dat$ctlR <- NULL

# Select motifs
ID <- readRDS("db/motif_counts/motifs_IDs.rds")
mot <- readRDS("db/motif_counts/twist008_motif_counts.rds")
mot <- melt(mot, "ID")

# Compute mean activity and residuals
motNeg <- mot[value==0, {
  dat[L %in% ID & R %in% ID, 
      .(act= mean(log2FoldChange),
        res= mean(residuals))]
}, variable]
motPos <- mot[value>=2, {
  dat[L %in% ID & R %in% ID, 
      .(act= mean(log2FoldChange),
        res= mean(residuals))]
}, variable]
res <- merge(motNeg, motPos, by= "variable", suffixes= c("_neg", "_pos"))
res[, `Mean activity [motif>1/motif=0] (log2)`:= act_pos-act_neg]
res[, `Mean residuals [motif>1/motif=0] (log2)`:= res_pos-res_neg]
res[ID, name:= i.cluster, on= "variable==motif_ID"]
res[, hit:= name %in% c("Trl/1", "Ebox/CATATG/twi", "DRE/1")]

pdf("pdf/draft/motif_impact_act_vs_res.pdf", 4, 4)
pl <- ggplot(res,
             aes(`Mean activity [motif>1/motif=0] (log2)`,
                 `Mean residuals [motif>1/motif=0] (log2)`,
                 color= hit,
                 label = name)) +
  scale_color_manual(values=c("grey50", "red"))+
  geom_point(size= 2,
             show.legend = FALSE) +
  geom_text_repel(max.overlaps= 10,
                  size= 4,
                  show.legend = FALSE) +
  theme(legend.position = "none")+
  theme_classic()
plot(pl)
dev.off()