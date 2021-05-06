setwd("/groups/stark/vloubiere/projects/pe_STARRSeq_2/")
# Packs ####
require(vlfunctions)
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
require(readxl)
require(PWMEnrich)
require(TFBSTools)
require(seqLogo)
require(motifmatchr)
require(kohonen)
####

# Processing ####
if(1==0)
{
  file.edit("git_peSTARRSeq/subscripts/extract_fastq_from_VBC_bam.R")
  file.edit("git_peSTARRSeq/subscripts/create_Rsubread_indexes.R")
  file.edit("git_peSTARRSeq/subscripts/alignment.R")
  file.edit("git_peSTARRSeq/subscripts/alignment_statistics.R")
  file.edit("git_peSTARRSeq/subscripts/compute_counts.R")
  file.edit("git_peSTARRSeq/subscripts/saturation_sequencing.R")
  file.edit("git_peSTARRSeq/subscripts/spike_ins_analyses.R")
  file.edit("git_peSTARRSeq/subscripts/template_swicthing_analyses.R")
  file.edit("git_peSTARRSeq/subscripts/umi_collapsing.R")
  file.edit("git_peSTARRSeq/subscripts/differential_analysis.R")
  file.edit("git_peSTARRSeq/subscripts/compute_expected_scores.R")
}

# Enhancer features ####
if(0==1)
{
  file.edit("git_peSTARRSeq/subscripts/compute_enhancer_features.R")
}

# Analyses ####
if(i==0)
{
  file.edit("git_peSTARRSeq/subscripts/cross_valid_bernardo.R")
  file.edit("git_peSTARRSeq/subscripts/modelling_activity.R")
  file.edit("git_peSTARRSeq/subscripts/modelling_residuals.R") # Not finished
  file.edit("git_peSTARRSeq/subscripts/close_vs_distant_pairs.R")
  file.edit("git_peSTARRSeq/subscripts/multiple_enhancer_loci_vs_few.R")
  file.edit("git_peSTARRSeq/subscripts/homotypic_vs_heterotypic_pairs.R")
  file.edit("git_peSTARRSeq/subscripts/heatmap_log2FC.R")
}

# Collaboration Bernardo ####
if(1==0)
{
  file.edit("git_peSTARRSeq/subscripts/clean_data_for_Bernardo.R")
}

# luciferase validations ####
if(1==0)
{
  file.edit("git_peSTARRSeq/subscripts/assemble_validations_metadata.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_analysis_luc_validations_constructs.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_barplot.R")
  file.edit("git_peSTARRSeq/subscripts/PCC_stripchart_luciferase_validations_peSTARR.R")
  file.edit("git_peSTARRSeq/subscripts/PCC_residuals_luciferase_validations_peSTARR.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_heatmap_pairs.R") # Not working yet
}










