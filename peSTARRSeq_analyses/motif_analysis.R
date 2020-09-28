setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(parallel)

dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[between(median_L, 1, 7) & between(median_R, 1, 7)]
dat <- dat[enh_L!=enh_R]

feat <- readRDS("Rdata/library/lib_features.rds")
feat <- melt(feat, id.vars = "ID", measure.vars = patterns("^motif__"))
feat <- feat[ID %in% dat$enh_L | ID %in% dat$enh_R]
feat[, variable:= gsub("motif__", "", variable)]
setkeyv(feat, c("ID", "variable"))

#-------------------------------------------#
# 1- Compute linear regression pairwise motif combinations
#-------------------------------------------#
if(!file.exists("Rdata/motifs/pairwise_motifs_combinations.rds"))
{
  combinations <- CJ(mot_L= unique(feat$variable), mot_R= unique(feat$variable))
  combinations <- combinations[,
                               {
                                 current <- data.table(feat[.(dat$enh_L, mot_L), .(enh_L= ID, counts_L= value), allow.cartesian= T],
                                                       feat[.(dat$enh_R, mot_R), .(enh_R= ID, counts_R= value), allow.cartesian= T])
                               }, .(mot_L, mot_R)]
  colc <- c("median_L", "median_R", "diff", "log2FC_add", "log2FoldChange") 
  combinations[dat, cols:= .(i.median_L, i.median_R, i.diff, i.log2FC_add, i.log2FoldChange), on= c("enh_L", "enh_R")]
  saveRDS(combinations, "Rdata/motifs/pairwise_motifs_combinations.rds")
}

if(!file.exists("Rdata/motifs/pairwise_motifs_combinations.rds"))
{
  combinations <- readRDS("Rdata/motifs/pairwise_motifs_combinations.rds")
  res <- combinations[,
                      {
                        res1 <- as.data.table(summary(lm(diff~median_L+median_R+counts_L+counts_R, .SD))$coefficients)
                        res2 <- as.data.table(summary(lm(diff~median_L+median_R+counts_L+counts_R+counts_L*counts_R, .SD))$coefficients)
                        res <- data.table(predictor= c("mot_L", "mot_R", "mot_L*mot_R"),
                                          tval= c(res1$`t value`[4:5], res2$`t value`[6]),
                                          pval= c(res1$`Pr(>|t|)`[4:5], res2$`Pr(>|t|)`[6]))
                      }, .(mot_L, mot_R)]
  saveRDS(res, "Rdata/motifs/lm_pairwise_motifs_combinations.rds")
}

if(!file.exists("Rdata/motifs/pairwise_motifs_combinations_LRIDs.rds"))
{
  combinations <- readRDS("Rdata/motifs/pairwise_motifs_combinations.rds")
  res <- combinations[,
                      {
                        res1 <- as.data.table(summary(lm(diff~enh_L+enh_R+counts_L+counts_R, .SD))$coefficients)
                        res2 <- as.data.table(summary(lm(diff~enh_L+enh_R+counts_L+counts_R+counts_L*counts_R, .SD))$coefficients)
                        res <- data.table(predictor= c("mot_L", "mot_R", "mot_L*mot_R"),
                                          tval= c(res1$`t value`[4:5], res2$`t value`[6]),
                                          pval= c(res1$`Pr(>|t|)`[4:5], res2$`Pr(>|t|)`[6]))
                      }, .(mot_L, mot_R)]
  saveRDS(res, "Rdata/motifs/lm_pairwise_motifs_combinations_LRIDs.rds")
}

#-------------------------------------------#
# 2- Heatmap
#-------------------------------------------#
# #import and simplify informative cluster names
motif_info <- readRDS("Rdata/motifs/som_enriched_motifs.rds")$info
motif_info[Dmel=="", Dmel:= best_match]
motif_info[grepl("^cisbp|^jaspar|^hocomoco|^predrem", Dmel), Dmel:= gsub("__", "", Dmel)]
motif_info <- motif_info[, .(name= unlist(tstrsplit(Dmel, "__"))), .(Dmel, best_match)]
motif_info[, is_seq:= gsub("A|C|G|T|R|Y|S|W|K|M|B|D|H|V|N", "", name)=="", name]
motif_info[, seq_only:= all(is_seq), Dmel]
motif_info[!(seq_only), name:= paste0(.SD[!(is_seq), name], collapse= ","), Dmel]
motif_info[(seq_only), name:= paste0(.SD[1, name], collapse= ","), Dmel]
motif_info <- unique(motif_info[, .(best_match, name)])
motif_info[, name:= gsub(" ", "_", name)]

res <- readRDS("Rdata/motifs/lm_pairwise_motifs_combinations.rds")
res[motif_info, Dmel_L:= i.name, on= "mot_L==best_match"]
res[motif_info, Dmel_R:= i.name, on= "mot_R==best_match"]
res[tval>0, pval:= -log10(pval)]
res[tval<0, pval:= log10(pval)]

t_L <- as.matrix(dcast(res[predictor=="mot_L"], Dmel_L~Dmel_R, value.var = "pval"), 1)
t_R <- as.matrix(dcast(res[predictor=="mot_R"], Dmel_L~Dmel_R, value.var = "pval"), 1)
t_LR <- as.matrix(dcast(res[predictor=="mot_L*mot_R"], Dmel_L~Dmel_R, value.var = "pval"), 1)

Cc <- colorRampPalette(c("cornflowerblue", "white", "white", "tomato"))(100)
br <- seq(-5, 5, length.out = 101)
t_LR[t_LR>5] <- 5
t_LR[t_LR< -5] <- -5
par(mar= c(18,18,5,7))
my_pheatmap(t_LR, color = Cc, breaks = br, row_labels = T, col_labels = T)




