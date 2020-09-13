setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/scripts/R_functions/R_shell_cmd_wrap.R")
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(GenomicRanges)
require(seqinr)
require(motifmatchr)
require(PWMEnrich) # module load gsl/2.1-foss-2017a
require(TFBSTools)
require(seqLogo)
require(digest)

lib <- readRDS("Rdata/C_features_final_table.rds")
lib <- lib[, uniq_ID:detail]

### Motifs counts  ############################
# load motifs and select the ones present in S2
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
DT <- as.data.table(TF_clusters_PWMs$metadata)
sel <- DT[S2_exp>1 | X..motif_collection_name %in% c("stark", "elemento"), motif_name]
sel <- match(sel, name(TF_clusters_PWMs$All_pwms_log_odds))

# Count hits low cutoff
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], GRanges(lib$coor), genome= "dm3", p.cutoff= 5e-4, bg="even", out= "scores")

counts <- as.matrix(motifCounts(hit))
colnames(counts) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_low_motif_count.", motif_name), on= "motif_name"]
rownames(counts) <- lib$uniq_ID
counts <- as.data.table(counts, keep.rownames = T)

scores <- as.matrix(motifmatchr::motifScores(hit))
colnames(scores) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_low_motif_score.", motif_name), on= "motif_name"]
rownames(scores) <- lib$uniq_ID
scores <- as.data.table(scores, keep.rownames = T)

# Add to lib
lib <- lib[counts, , on= "uniq_ID==rn"]
lib <- lib[scores, , on= "uniq_ID==rn"]

# Count hits high cutoff
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], GRanges(lib$coor), genome= "dm3", p.cutoff= 1e-4, bg="even", out= "scores")

counts <- as.matrix(motifCounts(hit))
colnames(counts) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_high_motif_count.", motif_name), on= "motif_name"]
rownames(counts) <- lib$uniq_ID
counts <- as.data.table(counts, keep.rownames = T)

scores <- as.matrix(motifmatchr::motifScores(hit))
colnames(scores) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_high_motif_score.", motif_name), on= "motif_name"]
rownames(scores) <- lib$uniq_ID
scores <- as.data.table(scores, keep.rownames = T)

# Add to lib
lib <- lib[counts, , on= "uniq_ID==rn"]
lib <- lib[scores, , on= "uniq_ID==rn"]
# melt
lib <- melt(lib, id.vars = colnames(lib[, uniq_ID:detail]), 
            measure.vars = patterns(low_motif_count= "_low_motif_count.", low_motif_score= "_low_motif_score.", 
                                    high_motif_count= "_high_motif_count.", high_motif_score= "_high_motif_score."))
lib[, motif:= name(TF_clusters_PWMs$All_pwms_log_odds[sel])[as.numeric(variable)]]
lib[as.data.table(TF_clusters_PWMs$metadata), i.Dmel, on= "motif==motif_name"]
lib <- merge(lib, as.data.table(TF_clusters_PWMs$metadata)[, .(motif_name, Dmel_prot= Dmel)], by.x= "motif", by.y= "motif_name", all.x= T)

########## Compute feat object
saveRDS(lib, "Rdata/D_motifs_final_table.rds")

### Work on subset of interest?
## Expression/prot cutoffs
# sel <- DT[S2_exp>1 & !is.na(S2_protein_exp) & !is.na(Dmel), motif_name]
# sel <- match(sel, name(TF_clusters_PWMs$All_pwms_log_odds))
## Curated for 1/ no redundancy 2/ enriched in enhancers OR known relevance
# sel <- c("homer__AVYTATCGATAD_DREF", "jaspar__MA0536.1", "cisbp__M2325", "idmmpmm__srp", "flyfactorsurvey__kay_Jra_SANGER_5_FBgn0001291", 
#          "flyfactorsurvey__Trl_FlyReg_FBgn0013263", "flyfactorsurvey__crp_SANGER_10_FBgn0001994", "cisbp__M4793", "stark__CACATGT", "cisbp__M6467", 
#          "homer__MYGGTCACACTG_Unknown1", "flyfactorsurvey__GATAd_SANGER_5_FBgn0032223", "flyfactorsurvey__HLH106_SANGER_10_FBgn0015234", "cisbp__M2575", 
#          "cisbp__M5212", "bergman__br-Z4", "stark__MAACAA", "flyfactorsurvey__jim_SANGER_2.5_FBgn0027339", "bergman__pnr")
# sel <- match(sel, name(TF_clusters_PWMs$All_pwms_log_odds))

## Alternative unbiased selection based on enrichment (misses AP-1 and few others...)
# counts <- melt(counts, "rn")
# counts <- merge(counts, lib[, .(rn= uniq_ID, twist_log2FC, twist_padj)])
# counts[, c("log2_OR", "fisher_pval") := 
#             {
#               current <- fisher.test(table(twist_log2FC>1 & twist_padj<0.00001, value>0))
#               list(log2(current$estimate), current$p.value)
#             }, variable]
# counts[, fisher_padj := p.adjust(fisher_pval)]
# counts <- counts[log2_OR>0 & fisher_pval<0.001]
## Count hits low cutoff
# hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], GRanges(lib$coor), genome= "dm3", p.cutoff= 1e-3, bg="even", out= "scores")
# counts <- as.matrix(motifCounts(hit))
# colnames(counts) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_lcutoff_motif.", motif_name), on= "motif_name"]
# rownames(counts) <- lib$uniq_ID
# counts <- as.data.table(counts, keep.rownames = T)
# # Add to lib
# lib <- lib[counts, , on= "uniq_ID==rn"]
# ## Count hits high cutoff
# hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], GRanges(lib$coor), genome= "dm3", p.cutoff= 1e-4, bg="even", out= "scores")
# counts <- as.matrix(motifCounts(hit))
# colnames(counts) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_hcutoff_motif.", motif_name), on= "motif_name"]
# rownames(counts) <- lib$uniq_ID
# counts <- as.data.table(counts, keep.rownames = T)
# # Add to lib
# lib <- lib[counts, , on= "uniq_ID==rn"]

