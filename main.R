setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

if(F)
{
  #-------------------------#
  # Libraries
  #-------------------------#
  # twist008 ---------------#
  file.edit("git_peSTARRSeq/subscripts/Generate_Ecoli_random_sequences.R")
  file.edit("git_peSTARRSeq/subscripts/OSC_specific_enhancers_bernardo.R")# Broken links
  file.edit("git_peSTARRSeq/subscripts/similarity_function.R") # used to select enh with diff seq twist008
  file.edit("git_peSTARRSeq/subscripts/design_twist_008.R") # lib can be found at "Rdata/vl_library_twist008_112019.rds"
  # twist012 ---------------# Original folder backed up at "old_versions/original_folder_twist12_design_backup/"
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_peak_calling.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_enrichment_quantif.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_candidates_selection.R")
  file.edit("git_peSTARRSeq/subscripts/CPs_selection_twist.R")
  file.edit("git_peSTARRSeq/subscripts/design_twist_012.R") # lib can be found at "Rdata/vl_library_twist12_210610.rds"
  
  #-------------------------#
  # Libraries features
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/make_300bp_uniq_enhancers_object.R")
  file.edit("git_peSTARRSeq/subscripts/motifs_clustering_300bp_enhancers.R")
  file.edit("git_peSTARRSeq/subscripts/final_lib_features.R")
  
  #-------------------------#
  # Luciferase assays
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/HAM_luciferase_assays.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_analyses.R")
  
  #-------------------------#
  # STARRSeq design
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/restriction_enzyme_selection.R") # NotI cuts the peSTARRSeq vector but is not detected!?
  file.edit("git_peSTARRSeq/subscripts/spacer_selection_revSTARRSeq.R")
  file.edit("git_peSTARRSeq/subscripts/primer_design_revSTARRSeq_spacers.R")
  file.edit("git_peSTARRSeq/subscripts/enh_promoter_distance.R")
  file.edit("git_peSTARRSeq/subscripts/qPCR_non_self_lig_vllib004.R")
  file.edit("git_peSTARRSeq/subscripts/intron_selection.R")
  file.edit("git_peSTARRSeq/subscripts/design_introns_primers.R")
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_to_replace_DSCP.R")
  
  #-------------------------#
  # Sanger sequencing
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/sanger_peSTARRSeq_SCR1_luc_validations.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_STARRSeq_backbone.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_actCP_STARRSeq_backbones.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_seqREADY_revSTARRSeq_single_col.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib004-007.R") # Check inverse PCR outcome, code broken
  
  #-------------------------#
  # peSTARR-Seq pipeline
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/pipeline_3.0.R")
  
  #-------------------------#
  # peSTARR-Seq modeling
  #-------------------------#
  
  

}
