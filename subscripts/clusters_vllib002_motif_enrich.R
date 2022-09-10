setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import Clustering object and retrieve enhancer sequences
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
cl <- readRDS("Rdata/clustering_lm_residuals_vllib002.rds")
rows <- cl$rows
rows[feat, enh_seq:= i.enh_seq, on= "name==ID"]
cols <- cl$cols
cols[feat, enh_seq:= i.enh_seq, on= "name==ID"]

# Select motifs using random forest?
# feat <- fread("Rdata/final_300bp_enhancer_features.txt")
# mot <- feat[rows$name, on= "ID"]
# mot <- mot[, which(names(mot) %in% vl_Dmel_motifs_DB_full$motif), with= F]
# set.seed(10)  # Setting seed
# rfL <- randomForest::randomForest(x = mot,
#                                   y = rows$cl,
#                                   ntree = 100,
#                                   importance= T)
# randomForest::varImpPlot(rfL,
#            sort = T,
#            n.var = 50,
#            main = "Top 10 - Variable Importance")


# Design random bg sequences
if(!file.exists("Rdata/clustering_lm_residuals_vllib002_motif_enrich.rds"))
{
  set.seed(1)
  rdm <- vl_random_regions_BSgenome(genome = "dm3", 
                                    n= 1000,
                                    width = 249)
  rdm <- BSgenome::getSeq(BSgenome::getBSgenome("dm3"), 
                          rdm$seqnames, 
                          rdm$start, 
                          rdm$end, 
                          as.character= T)
  
  # Compute motifs
  res <- list()
  rdm_mot <- vl_motif_counts(sequences = rdm)
  for(cdition in c("rows", "cols"))
  {
    mot <- rbind(vl_motif_counts(sequences = get(cdition)$enh_seq),
                 rdm_mot)
    res[[cdition]] <- vl_motif_cl_enrich(mot,
                                         cl_IDs = c(as.character(get(cdition)$cl), rep("rdm", nrow(rdm_mot))),
                                         collapse_clusters = unlist(tstrsplit(colnames(mot), "__", keep= 1)), 
                                         control_cl = "rdm", 
                                         plot= F)
  }
  names(res) <- c("L", "R")
  saveRDS(res, 
          "Rdata/clustering_lm_residuals_vllib002_motif_enrich.rds")
}

##################################################
# PLOT
##################################################
res <- readRDS("Rdata/clustering_lm_residuals_vllib002_motif_enrich.rds")
res$L[, cl:= factor(cl, c("No syn.", "Weak syn.", "Medium syn.", "High syn."))]
res$R[, cl:= factor(cl, c("No syn.", "Weak syn.", "Medium syn.", "High syn."))]


.c <- na.omit(res$L)
setorderv(.c, "padj")
sel <- .c[, variable[log2OR>0][1:15], cl]$V1
dmat <- dcast(.c[variable %in% sel], 
              variable~cl, 
              value.var = "log2OR")
dmat <- as.matrix(dmat, 1)
dmat <- apply(dmat, 2, function(x) {
  x[x==Inf] <- max(x[x!=Inf])
  x[x==(-Inf)] <- min(x[x!=(-Inf)])
  return(x)
})
par(mgp= c(1.5, 0.5, 0),
    mar= c(5,8,1,6),
    tcl= -0.2,
    las= 1,
    cex= 0.6)
hm <- vl_heatmap(dmat, 
                 cluster_cols= F)
grep("DRE", rownames(dmat))


pdf("pdf/draft/cluster_motifs_enrichment_vllib002.pdf", 
    height = 4.5, 
    width = 4)
par(mgp= c(1.5, 0.5, 0),
    mar= c(2,10,1,6),
    tcl= -0.2,
    las= 1,
    cex= 0.6)
plot(res$L,
     padj_cutoff= 0.01,
     cex.balloons= 0.5, 
     top_enrich= 20,
     col= c("blue", "red"))
plot(res$R,
     padj_cutoff= 0.01,
     cex.balloons= 0.5, 
     top_enrich= 20, 
     col= c("blue", "red"))
dev.off()
