setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Select motifs
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
sel <- as.data.table(TF_clusters_PWMs[["metadata"]])[S2_exp>0, motif_name]

#--------------------------------------------#
# Random control set
#--------------------------------------------#
set.seed(1)
ctl <- vl_random_regions_BSgenome(genome = "dm3", 
                                  n= 1000,
                                  width = 249)
# Low stringency counts
counts <- vl_motif_counts(bed = ctl, 
                          sel = sel,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
counts[, ID:= ctl[, paste0(seqnames, ":", start, "-", end, ":", strand)]]
setcolorder(counts, "ID")
saveRDS(counts, 
        "db/motif_counts/random_controls_1000_low_stringency_no_collapsing.rds")

#--------------------------------------------#
# twist 008
#--------------------------------------------#
# Low stringency counts
lib <- readRDS("Rdata/vl_library_twist008_112019.rds")
counts <- vl_motif_counts(sequences = lib$enh_sequence, 
                          sel = sel,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
counts[, ID:= lib$ID_vl]
setcolorder(counts, "ID")
saveRDS(counts, 
        "db/motif_counts/twist008_motif_counts_low_stringency_no_collapsing.rds")

#--------------------------------------------#
# twist 012
#--------------------------------------------#
# Low stringency counts
lib <- readRDS("Rdata/vl_library_twist12_210610.rds")
counts <- vl_motif_counts(sequences = lib$enh_seq, 
                          sel = sel,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
counts[, ID:= lib$ID]
setcolorder(counts, "ID")
saveRDS(counts, 
        "db/motif_counts/twist012_motif_counts_low_stringency_no_collapsing.rds")
