setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Select motifs
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
sel <- as.data.table(TF_clusters_PWMs[["metadata"]])
sel[, check:= any(S2_exp>=5), Motif_cluster_name]
sel <- sel[(check), motif]

# Indeitify relevant motifs ----
## Control set ----
set.seed(1)
rdm <- vl_random_regions_BSgenome(genome = "dm3", n= 1000, width = 249)
ctl <- vl_getSequence(rdm, "dm3")
lib8 <- readRDS("Rdata/vl_library_twist008_112019.rds")$enh_sequence
lib12 <- readRDS("Rdata/vl_library_twist12_210610.rds")$enh_seq
lib15 <- readRDS("Rdata/vl_library_twist015_112022.rds")$enh_sequence
## Counts ----
dat <- list(ctl= ctl,
            lib8= lib8,
            lib12= lib12,
            lib15= lib8)
dat <- lapply(dat, as.data.table)
dat <- rbindlist(dat, idcol = "lib")
counts <- vl_motif_counts(sequences = dat$V1, 
                          sel = sel,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
enr <- vl_motif_cl_enrich(split(counts, dat$lib),
                          control_cl = "ctl",
                          plot = F)
enr[vl_Dmel_motifs_DB_full, cluster:= i.motif_cluster, on= "variable==motif_ID"]
enr[vl_Dmel_motifs_DB_full, Dmel:= i.Dmel, on= "variable==motif_ID"]
sel <- enr[padj<1e-5 & log2OR>0 & set_hit>10, .SD[which.max(log2OR), .(Dmel, motif_ID= variable)], cluster]
saveRDS(sel, "db/motif_counts/motifs_IDs.rds")

# Random controls ----
counts <- vl_motif_counts(bed = rdm, 
                          sel = sel$motif_ID,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
# setnames(counts, new= as.character(sel$cluster))
counts[, ID:= rdm[, paste0(seqnames, ":", start, "-", end, ":", strand)]]
setcolorder(counts, "ID")
saveRDS(counts, 
        "db/motif_counts/random_controls_1000.rds")

# Twist 008 ----
lib <- readRDS("Rdata/vl_library_twist008_112019.rds")
counts <- vl_motif_counts(sequences = lib$enh_sequence, 
                          sel = sel$motif_ID,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
# setnames(counts, new= as.character(sel$cluster))
counts[, ID:= lib$ID_vl]
setcolorder(counts, "ID")
saveRDS(counts, 
        "db/motif_counts/twist008_motif_counts.rds")

# twist 012 ----
lib <- readRDS("Rdata/vl_library_twist12_210610.rds")
counts <- vl_motif_counts(sequences = lib$enh_seq, 
                          sel = sel$motif_ID,
                          genome= "dm3",
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
# setnames(counts, new= as.character(sel$cluster))
counts[, ID:= lib$ID]
setcolorder(counts, "ID")
saveRDS(counts, 
        "db/motif_counts/twist012_motif_counts.rds")

# twist 015 ----
lib <- readRDS("Rdata/vl_library_twist015_112022.rds")
counts <- vl_motif_counts(sequences = lib$enh_seq, 
                          sel = sel$motif_ID,
                          genome= BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3,
                          bg = "genome", 
                          p.cutoff = 5e-04, 
                          collapse_overlapping = F)
# setnames(counts, new= as.character(sel$cluster))
counts[, ID:= lib$ID]
setcolorder(counts, "ID")
saveRDS(counts,
        "db/motif_counts/twist015_motif_counts.rds")
