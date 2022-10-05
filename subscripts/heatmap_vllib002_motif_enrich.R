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
feat[cl$rows, clL:= i.cl, on= "ID==name"]
feat[cl$cols, clR:= i.cl, on= "ID==name"]
feat[, cl:= paste0(clL, "__", clR)]
feat <- feat[!is.na(clL) & !is.na(clR)]

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


res <- vl_motif_cl_enrich(rbind(as.matrix((feat[, top$rn, with= F])), 
                                motCtl),
                           cl_IDs= factor(c(as.character(feat$cl), 
                                            rep("rdm", nrow(motCtl)))),
                           collapse_clusters = vl_Dmel_motifs_DB_full[colnames(motL), motif_cluster, on= "motif"],
                           control_cl = "rdm")

test <- res[variable=="DRE/1" & cl!= "rdm"]
test[, c("clL", "clR"):= tstrsplit(cl, "__")]
test[, padj:= ifelse(log2OR>0, -log10(padj), log10(padj))]
test <- dcast(test, clL~clR, value.var = "log2OR")
test <- as.matrix(test, 1)
test[test==(-Inf)] <-0
vl_heatmap(test, cluster_rows = F, cluster_cols = F)