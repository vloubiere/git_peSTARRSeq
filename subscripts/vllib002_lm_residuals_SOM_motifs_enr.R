setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import Clustering object and retrieve enhacner sequences
cl <- readRDS("Rdata/vllib002_lm_residuals_SOM.rds")
rows <- cl$rows
cols <- cl$cols
feat <- fread("Rdata/final_300bp_enhancer_features.txt")[, .(ID, enh_seq)]
rows[feat, enh_seq:= i.enh_seq, on= "name==ID"]
cols[feat, enh_seq:= i.enh_seq, on= "name==ID"]

# Design random bg sequences
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
        "Rdata/vllib002_lm_residuals_SOM_motifs_enr.rds")