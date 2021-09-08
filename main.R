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
  file.edit("git_peSTARRSeq/subscripts/informative_motif_counts.R") # I do not use som clustering anymore (top 50 motifs)
  file.edit("git_peSTARRSeq/subscripts/SOM_informative_motif_counts.R") # SOM
  file.edit("git_peSTARRSeq/subscripts/chromatin_features_and_gene_assignment.R")
  file.edit("git_peSTARRSeq/subscripts/final_lib_features.R")
  
  #-------------------------#
  # Luciferase assays
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/HAM_luciferase_assays.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_vs_STARRSeq_residuals.R")
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_PCC_and_connected_stripchart.R")
  
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
  file.edit("git_peSTARRSeq/subscripts/sanger_HAM_SUP1_luciferase_constructs.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_peSTARRSeq_SCR1_luc_validations.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_STARRSeq_backbone.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_actCP_STARRSeq_backbones.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_seqREADY_revSTARRSeq_single_col.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib004-007.R") # Check inverse PCR outcome, code broken
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib015-016.R") # First libs dumped -> no gel purif before gibson, many single enhancers
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib015-016_#2.R") # Second libs OK
  
  #-------------------------#
  # peSTARR-Seq pipeline
  #-------------------------#
  file.edit("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/create_twist8_subread_index.R")
  file.edit("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/create_twist12_subread_index.R")
  file.edit("git_peSTARRSeq/subscripts/pipeline_3.0.R")
  file.edit("git_peSTARRSeq/subscripts/basic_alignment_statistics.R") # Used by the pipeline
  file.edit("git_peSTARRSeq/subscripts/aggregate_alignment_statistics.R") # Used by the pipeline
  file.edit("git_peSTARRSeq/subscripts/sequencing_saturation.R")
  file.edit("git_peSTARRSeq/subscripts/PCC.R")
  file.edit("git_peSTARRSeq/functions/generate_final_data_table.R")
  
  #-------------------------#
  # Aggregate activity
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/heatmap_activity_aggreagte.R")
  file.edit("git_peSTARRSeq/subscripts/smoothsScatter_obs_exp.R")
  file.edit("git_peSTARRSeq/subscripts/boxplot_FC_silencer_insulator_LR.R")
  
  #-------------------------#
  # Fine grained heatmaps
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/heatmap_activity.R")
  file.edit("git_peSTARRSeq/subscripts/detailed_heatmap_additivity.R")
  
  #-------------------------#
  # peSTARR-Seq modeling
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/linear_models_activity_residuals.R")
  file.edit("git_peSTARRSeq/subscripts/smoothScatterplot_motifs_contribution.R")
  file.edit("git_peSTARRSeq/subscripts/barplot_motifs_contribution.R")
  
  #-------------------------#
  # Spacer length impact
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/comparison_vllib002_vllib006.R")
  
  #-------------------------#
  # Others
  #-------------------------#
  # shn screenshot
  file.edit("git_peSTARRSeq/subscripts/shn_screenshot.R")
  # Average CHIP-Seq tracks
  file.edit("git_peSTARRSeq/subscripts/aggregate_ChIPSeq_tracks.R")
  # Get 200bp STARR-Seq peaks from Bernardo
  file.edit("git_peSTARRSeq/subscripts/get_200bp_STARRSeq_peaks_BA.R")
  # BA clusters motif enrichment
  file.edit("git_peSTARRSeq/subscripts/get_clusters_motifs_enrichment_BA.R")

  #-------------------------#
  # OLD
  #-------------------------#
}
