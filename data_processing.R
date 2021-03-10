setwd("/groups/stark/vloubiere/projects/pe_STARRSeq_2/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
# Packs ####
require(data.table)
require(Rsubread)
require(Biostrings)
require(parallel)
require(gridExtra)
require(pheatmap)
require(DESeq2)
require(seqinr)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(motifmatchr)
require(PWMEnrich)
require(TFBSTools)
require(seqLogo)
require(motifmatchr)
require(kohonen)
####

# Source ####
if(1==0)
{
  source("git_peSTARRSeq/subscripts/extract_fastq_from_VBC_bam.R")
  source("git_peSTARRSeq/subscripts/create_Rsubread_indexes.R")
  source("git_peSTARRSeq/subscripts/compute_enhancer_features.R")
  source("git_peSTARRSeq/subscripts/alignment.R")
  source("git_peSTARRSeq/subscripts/alignment_statistics.R")
  source("git_peSTARRSeq/subscripts/compute_counts.R")
  source("git_peSTARRSeq/subscripts/saturation_sequencing.R")
  source("git_peSTARRSeq/subscripts/spike_ins_analyses.R")
  source("git_peSTARRSeq/subscripts/template_swicthing_analyses.R")
  source("git_peSTARRSeq/subscripts/differential_analysis.R")
  source("git_peSTARRSeq/subscripts/compute_expected_scores.R")
  source("git_peSTARRSeq/subscripts/cross_valid_bernardo.R")
  source("git_peSTARRSeq/subscripts/modelling_activity.R")
  source("git_peSTARRSeq/subscripts/modelling_residuals.R")
  source("git_peSTARRSeq/subscripts/close_vs_distant_pairs.R")
  source("git_peSTARRSeq/subscripts/multiple_enhancer_loci_vs_few.R")
  source("git_peSTARRSeq/subscripts/homotypic_vs_heterotypic_pairs.R")
  source("git_peSTARRSeq/subscripts/heatmap_log2FC.R")
  
  source("git_peSTARRSeq/subscripts/clean_data_for_Bernardo.R")
}





