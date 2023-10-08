setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import library sequences ---
lib8 <- readRDS("Rdata/vl_library_twist008_112019.rds")
lib8 <- as.data.table(lib8)
setnames(lib8,
         c("ID_vl", "enh_sequence"),
         c("ID", "sequence"))

# Random control set ----
set.seed(1)
rdm <- vl_random_regions_BSgenome(genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3, n= 1000, width = 249)
ctl <- vl_getSequence(rdm, BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3)
ctl <- data.table(ID= paste0("rdm_", rdm[, paste0(seqnames, ":", start, "-", end, ":", strand)]),
                  sequence= ctl,
                  group= "rdm")

# Merge the two groups ----
dat <- rbind(lib8[, .(group, ID, sequence)],
             ctl)

# Select motifs whose cognate TFs are expressed in S2 cells
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
sel <- as.data.table(TF_clusters_PWMs[["metadata"]])
sel[, check:= any(S2_exp>=5), Motif_cluster_name]
sel <- sel[(check), motif_name]

# Indentify motifs enriched in the oligo pool compared to ramdom sequences ----
counts <- vl_motif_counts(sequences = dat$sequence, 
                          sel = sel,
                          genome= BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3,
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
enr <- vl_motif_cl_enrich(split(counts, dat$group),
                          control_cl = "rdm")
enr[vl_Dmel_motifs_DB_full, cluster:= i.motif_cluster, on= "variable==motif_ID"]
enr[vl_Dmel_motifs_DB_full, Dmel:= i.Dmel, on= "variable==motif_ID"]
sel <- enr[padj<1e-5 & log2OR>1 & set_hit>20, .SD[which.max(log2OR), .(Dmel, motif_ID= variable)], cluster]
saveRDS(sel,
        "db/motif_counts/lib8_motifs_IDs.rds")

# Save counts motifs of interest ----
counts <- counts[, match(sel$motif_ID, names(counts)), with= F]
rdm <- counts[dat$group=="rdm"]
rdm[, ID:= dat[group=="rdm", ID]]
setcolorder(rdm, "ID")
saveRDS(rdm, 
        "db/motif_counts/random_controls_1000.rds")
enh <- counts[dat$group!="rdm"]
enh[, ID:= dat[group!="rdm", ID]]
setcolorder(enh, "ID")
saveRDS(enh, 
        "db/motif_counts/twist008_motif_counts.rds")
