setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Retrieve sequences ----
if(!file.exists("db/motif_counts/lib8_sequences_and_random_controls.rds"))
{
  ## Import library sequences ----
  lib8 <- readRDS("Rdata/vl_library_twist008_112019.rds")
  lib8 <- as.data.table(lib8)
  setnames(lib8,
           c("ID_vl", "enh_sequence"),
           c("ID", "sequence"))
  
  ## Random control set ----
  set.seed(1)
  rdm <- vl_random_regions_BSgenome(genome = BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3, n= 1000, width = 249)
  ctl <- vl_getSequence(rdm, BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3)
  ctl <- data.table(ID= paste0("rdm_", rdm[, paste0(seqnames, ":", start, "-", end, ":", strand)]),
                    sequence= ctl,
                    group= "rdm")
  
  ## Merge the two groups ----
  dat <- rbind(lib8[, .(group, ID, sequence)],
               ctl)
  saveRDS(dat,
          "db/motif_counts/lib8_sequences_and_random_controls.rds")
}else
  dat <- readRDS("db/motif_counts/lib8_sequences_and_random_controls.rds")

# Motif counts ----
if(!file.exists("db/motif_counts/twist008_motif_counts_full_table.rds"))
{
  ## Import motifs ----
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  sel <- as.data.table(TF_clusters_PWMs[["metadata"]])
  
  ## Select motifs whose cognate TFs are expressed in S2/OSC/ECD cells ----
  rep1 <- fread("db/public/GSM1000414_OSC_RNAseq_rep1.rpkm.txt.gz")
  rep2 <- fread("db/public/GSM1000415_OSC_RNAseq_rep2.rpkm.txt.gz")
  OSC <- merge(rep1, rep2, by= "V1")[V4.x>=5 & V4.y>=5, V1]
  rep1 <- fread("db/public/GSM1357052_S2_RNAseq_with_ecd_rep1.rpkm.txt.gz")
  rep2 <- fread("db/public/GSM1357053_S2_RNAseq_with_ecd_rep2.rpkm.txt.gz")
  ECD <- merge(rep1, rep2, by= "V1")[V4.x>=5 & V4.y>=5, V1]
  sel[, check:= any(S2_exp>=5 | FBgn %in% OSC | FBgn %in% ECD), Motif_cluster_name]
  sel <- sel[(check), motif_name]
  
  ## Motif counts
  counts <- vl_motif_counts(sequences = dat$sequence, 
                            sel = sel,
                            genome= BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3,
                            bg = "genome", 
                            p.cutoff = 5e-04)
  counts[, ID:= dat$ID]
  setcolorder(counts, "ID")
  saveRDS(counts,
          "db/motif_counts/twist008_motif_counts_full_table.rds")
  counts$ID <- NULL
}else
  counts <- readRDS("db/motif_counts/twist008_motif_counts_full_table.rds")[, -1]

# Identify motifs enriched in the oligo pool compared to ramdom sequences ----
if(!file.exists("db/motif_counts/lib8_motifs_enrichments.rds"))
{
  ## Enrichment ----
  enr <- vl_motif_cl_enrich(split(counts, dat$group),
                            control.cl = "rdm")
  setnames(enr,
           c("variable", "name"),
           c("motif_ID", "cluster"))
  ## Make clean names  for motifs of interest ----
  enr[, name:= cluster]
  enr[name=="GATA/1", name:= "GATA"]
  enr[name=="AP1/1", name:= "AP-1"]
  enr[name=="SREBP/1", name:= "SREBP"]
  enr[name=="DRE/1", name:= "Dref.1"]
  enr[name=="DRE/2", name:= "Dref.2"]
  enr[name=="Ebox/CATATG/twi", name:= "Twist"]
  enr[name=="Trl/1", name:= "Trl"]
  enr[name=="NR/7", name:= "EcR/usp/Hr3"]
  enr[name=="MAF/2", name:= "MAF_tj"]
  # Add Pwms ----
  enr[vl_Dmel_motifs_DB_full, pwm:= i.pwms_perc, on="motif_ID"]
  enr[, pwm:= lapply(pwm, as.matrix)]
  ## Save ----
  saveRDS(enr,
          "db/motif_counts/lib8_motifs_enrichments.rds")
}else
  enr <- readRDS("db/motif_counts/lib8_motifs_enrichments.rds")

# Save slected motifs ----
selMot <- enr[padj<1e-5 & set_hit>=5, .SD[which.min(padj)], .(cl, cluster)]
selMot <- selMot[, .SD[1], motif_ID]
selMot$cl <- selMot$pval <- selMot$log2OR <- selMot$padj <- NULL

selCounts <- counts[dat$group!="rdm", selMot$motif_ID, with= F]
selCounts[, ID:= dat[group!="rdm", ID]]
setcolorder(selCounts, "ID")
saveRDS(selCounts, 
        "db/motif_counts/twist008_motif_counts_selected.rds")

selCounts <- counts[dat$group=="rdm", selMot$motif_ID, with= F]
selCounts[, ID:= dat[group=="rdm", ID]]
setcolorder(selCounts, "ID")
saveRDS(selCounts, 
        "db/motif_counts/twist008_motif_counts_selected_negative_controls.rds")