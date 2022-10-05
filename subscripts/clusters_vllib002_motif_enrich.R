setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Top lasso motifs
models <- readRDS("Rdata/LASSO_residuals_vllib002.rds")
top <- rbind(models$L$coeffs, models$R$coeffs)
top <- top[vl_Dmel_motifs_DB_full, motif_cluster:= motif_cluster, on= "rn==motif"]
top <- na.omit(top)[order(abs(s0), decreasing = T)][s0!=0]
top <- top[, .SD[1], motif_cluster]

# Motif matrices
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
cl <- readRDS("Rdata/clustering_lm_residuals_vllib002.rds")
motL <- feat[cl$rows, top$rn, on= "ID==name", with= F]
motR <- feat[cl$cols, top$rn, on= "ID==name", with= F]

# Control matrix
set.seed(1)
rdm <- vl_random_regions_BSgenome(genome = "dm3", 
                                  n= 1000,
                                  width = 249)
rdm <- BSgenome::getSeq(BSgenome::getBSgenome("dm3"), 
                        rdm$seqnames, 
                        rdm$start, 
                        rdm$end, 
                        as.character= T)
motCtl <- vl_motif_counts(sequences = rdm, sel= top$rn)

resL <- vl_motif_cl_enrich(rbind(motL, motCtl),
                           cl_IDs= factor(c(as.character(cl$rows$cl), 
                                            rep("rdm", nrow(motCtl))),
                                          levels= c("rdm", "High syn.", "Medium syn.", "Weak syn.", "No syn.")),
                           collapse_clusters = vl_Dmel_motifs_DB_full[colnames(motL), motif_cluster, on= "motif"],
                           control_cl = "rdm",
                           plot= F)
resR <- vl_motif_cl_enrich(rbind(motR, motCtl),
                           cl_IDs= factor(c(as.character(cl$cols$cl), 
                                            rep("rdm", nrow(motCtl))),
                                          levels= c("rdm", "High syn.", "Medium syn.", "Weak syn.", "No syn.")),
                           collapse_clusters = vl_Dmel_motifs_DB_full[colnames(motR), motif_cluster, on= "motif"],
                           control_cl = "rdm",
                           plot= F)

######################################################
# PLOT
######################################################
pdf("pdf/draft/cluster_motifs_enrichment_vllib002.pdf", 
    height = 4.5, 
    width = 4.2)
par(mgp= c(1.5, 0.5, 0),
    mar= c(5.5,25,2,6),
    tcl= -0.2,
    las= 2,
    cex= 0.6)
plL <- plot(resL, 
            padj_cutoff= 5e-2,
            col= c("blue", "red"), 
            top_enrich= 10,
            cex.balloons= 0.6,
            main = "5' clusters")
vl_add_motifs(DT = plL$DT, cex.height = 1.2)
plR <- plot(resR, 
     padj_cutoff= 5e-2,
     col= c("blue", "red"), 
     top_enrich= 10,
     cex.balloons= 0.6,
     main = "3' clusters")
vl_add_motifs(plR$DT, cex.height = 1.2)
dev.off()
