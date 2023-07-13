setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

#######################################################################################################
# pe-STARR-Seq PIPELINE 
#######################################################################################################
if(F)
{
  # Pipeline 3.0 ------------------------------------------------------------------------#
  file.edit("git_peSTARRSeq/subscripts/create_twist8_subread_index.R")
  file.edit("git_peSTARRSeq/subscripts/create_twist12_subread_index.R")
  file.edit("git_peSTARRSeq/subscripts/create_twist15_subread_index.R")
  
  file.edit("git_peSTARRSeq/subscripts/pipeline_peSTARRSeq.R")
  file.edit("git_peSTARRSeq/subscripts/run_pipeline_parallel.R")
  
  # QC and sanity check -----------------------------------------------------------------#
  file.edit("git_peSTARRSeq/subscripts/Benchmarking_UMI_collapsing.R") # Compare to BA approach
  file.edit("git_peSTARRSeq/subscripts/Compare_DESeq2_to_raw_log2_ratio.R") # Sanity check
  file.edit("git_peSTARRSeq/subscripts/density_reads_position.R") # Sanity check
  file.edit("git_peSTARRSeq/subscripts/aggregate_alignment_statistics.R")
  file.edit("git_peSTARRSeq/subscripts/sequencing_saturation.R")
  file.edit("git_peSTARRSeq/subscripts/barplot_read_per_theoretical_pair.R")
  file.edit("git_peSTARRSeq/subscripts/PCC.R") 
}

#######################################################################################################
# PAPER ANALYSES 
#######################################################################################################
# Wrappers to clean folder
file.edit("git_peSTARRSeq/subscripts/clean_scripts_pdf.R")
if(F)
{
  # Objects
  file.edit("git_peSTARRSeq/subscripts/Libraries_motif_counts.R")
  file.edit("git_peSTARRSeq/subscripts/linear_model_vllib002.R")
  file.edit("git_peSTARRSeq/subscripts/PCA_vllib002_actPairs_residuals.R")
  file.edit("git_peSTARRSeq/subscripts/lasso_model_vllib002.R")
  file.edit("git_peSTARRSeq/subscripts/motif_enrichment_residuals.R")
  file.edit("git_peSTARRSeq/subscripts/Make_final_tables.R")
  
  # Alignment statistics
  file.edit("git_peSTARRSeq/subscripts/alignment_statistics.R")
  # PCC between replicates
  file.edit("git_peSTARRSeq/subscripts/PCC_peSTARRSeq_replicates.R")
  # PCC pSTARR-Seq/luciferase
  file.edit("git_peSTARRSeq/subscripts/Luciferase_validations.R")
  # Left vs right activities
  file.edit("git_peSTARRSeq/subscripts/Correlation_left_rigth_activities.R")
  # Single enhancers vs pairs
  file.edit("git_peSTARRSeq/subscripts/Candidates_classification_vllib002.R")
  file.edit("git_peSTARRSeq/subscripts/Compare_individual_vs_enh_pairs.R")
  # Activity heatmap ordered by activity
  file.edit("git_peSTARRSeq/subscripts/heatmap_act_ordered_ind_act_dev_pairs_vllib002.R")
  # Modelling vllib002 
  file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_vllib002.R")
  file.edit("git_peSTARRSeq/subscripts/smoothScatter_vllib002_actPairs_lm.R")
  # CP saturation by strong enhancers
  file.edit("git_peSTARRSeq/subscripts/promoter_saturation.R")
  # Modelling residuals 
  file.edit("git_peSTARRSeq/subscripts/evaluate_lasso_model.R")
  # Residuals heatmap ordered by residuals
  file.edit("git_peSTARRSeq/subscripts/heatmap_residuals_ordered_pca_vllib002.R")
  
  file.edit("git_peSTARRSeq/subscripts/linear_model_vllib015_vllib016.R")
  
  #------------------------------------#
  # Not used
  #------------------------------------#
  # Saturation
  file.edit("git_peSTARRSeq/subscripts/umi_saturation.R") # No saturation of the UMI
  file.edit("git_peSTARRSeq/subscripts/dev_enhancer_stength_vs_cooperativity.R") # Strong enhancers dont synergize and dont show specific motifs
  # Analysis vllib002 residuals Trl/Twist
  file.edit("git_peSTARRSeq/subscripts/lm_residuals_examples_dev_pairs_vllib002.R") # Linear model examples
  file.edit("git_peSTARRSeq/subscripts/density_residuals_motif_enrichment_dev_pairs_vllib002.R") # Show activity/residuals relationship
  file.edit("git_peSTARRSeq/subscripts/heatmaps_motifs_residuals_dev_pairs_vllib002.R") # Heatmap lm coeffs to predict residuals
  file.edit("git_peSTARRSeq/subscripts/activity_predictive_motifs_activity_matched_controls.R") # Boxplot ind vs comb for act matched pairs
  # Mutant library
  file.edit("git_peSTARRSeq/subscripts/mutant_library.R")
  file.edit("git_peSTARRSeq/subscripts/mutant_library_sequencing_depth.R")
  # DHS library
  file.edit("git_peSTARRSeq/subscripts/linear_model_vllib030.R") # train linear model on the whole library
  file.edit("git_peSTARRSeq/subscripts/DHS_library.R")
}

#######################################################################################################
# LIBRARY DESIGN AND QCs
#######################################################################################################
if(F)
{
  # Libraries ---------------------------------------------------------------------------#
  # twist008 
  file.edit("git_peSTARRSeq/subscripts/Generate_Ecoli_random_sequences.R")
  file.edit("git_peSTARRSeq/subscripts/OSC_specific_enhancers_bernardo.R")# Broken links
  file.edit("git_peSTARRSeq/subscripts/similarity_function.R") # used to select enh with diff seq twist008
  file.edit("git_peSTARRSeq/subscripts/design_twist_008.R") # lib can be found at "Rdata/vl_library_twist008_112019.rds"
  # twist012 -> Original folder backed up at "old_versions/original_folder_twist12_design_backup/"
  file.edit("git_peSTARRSeq/subscripts/SUHW_peak_selection.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_peak_calling.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_enrichment_quantif.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_candidates_selection.R")# Somehow, peaks are different from original run. See old folders
  file.edit("git_peSTARRSeq/subscripts/CPs_selection_twist.R")
  file.edit("git_peSTARRSeq/subscripts/design_twist_012.R") # lib can be found at "Rdata/vl_library_twist12_210610.rds"
  file.edit("git_peSTARRSeq/subscripts/aggregate_ChIPSeq_tracks_SUHW_silencers.R")
  # twist015
  file.edit("git_peSTARRSeq/subscripts/design_mutant_sequences_twist_015.R")
  file.edit("git_peSTARRSeq/subscripts/predict_mutant_sequences_twist_015.sh") #Not working cause I would need Bernies python environment
  file.edit("git_peSTARRSeq/subscripts/select_mutant_sequences_twist_015.R")
  file.edit("git_peSTARRSeq/subscripts/select_DHS_sites_twist_015.R")
  file.edit("git_peSTARRSeq/subscripts/design_twist_015.R")
  
  # STARRSeq reporter -------------------------------------------------------------------#
  # revSTARRSeq
  file.edit("git_peSTARRSeq/subscripts/restriction_enzyme_selection.R") # NotI cuts the peSTARRSeq vector but is not detected!?
  file.edit("git_peSTARRSeq/subscripts/spacer_selection_revSTARRSeq.R")
  file.edit("git_peSTARRSeq/subscripts/primer_design_revSTARRSeq_spacers.R")
  file.edit("git_peSTARRSeq/subscripts/enh_promoter_distance.R")
  file.edit("git_peSTARRSeq/subscripts/qPCR_non_self_lig_vllib004.R")
  file.edit("git_peSTARRSeq/subscripts/compare_vllib006_protocols.R")
  # Introns 
  file.edit("git_peSTARRSeq/subscripts/intron_selection.R") # without enhancer
  file.edit("git_peSTARRSeq/subscripts/intron_enhancer_selection.R") # with enhancer
  file.edit("git_peSTARRSeq/subscripts/design_introns_enhancer_primers.R")
  # Change promoter
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_to_replace_DSCP.R") # Mistake! Some CPs contain AgeI/SalI sites (see next line)
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_no_enzymatic_restriciton_to_replace_DSCP.R") # Not optimal cause I overlooked hk induction!
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_to_replace_DSCP_final.R")
  
  # Sanger sequencing -------------------------------------------------------------------#
  # Backbones
  file.edit("git_peSTARRSeq/subscripts/sanger_actCP_STARRSeq_backbones.R") # Failed, the CP is DSCP
  file.edit("git_peSTARRSeq/subscripts/sanger_actCP3_STARRSeq_backbone.R") # Failed, the CP is DSCP
  file.edit("git_peSTARRSeq/subscripts/sanger_alternative_CPs_STARRSeq_final.R")
  # Libraries single cols
  file.edit("git_peSTARRSeq/subscripts/sanger_seqREADY_revSTARRSeq_single_col.R")
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib004-007.R") # Check inverse PCR outcome, code broken
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib015-016.R") # First libs dumped -> no gel purif before gibson, many single enhancers
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib015-016_#2.R") # Second libs OK
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib017-020.R") # Seems OK
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib021-022.R") # Seems OK
  file.edit("git_peSTARRSeq/subscripts/sanger_vllib023-028.R") # Seems OK
}

#######################################################################################################
# Luciferase validations 
#######################################################################################################
if(F)
{
  # HAM CP (Lorena) ---------------------------------------------------------------------#
  # Sanger seq
  file.edit("git_peSTARRSeq/subscripts/sanger_HAM_SUP1_luciferase_constructs.R")
  # Luciferase
  file.edit("git_peSTARRSeq/subscripts/HAM_luciferase_assays.R") # Luciferase
  
  # dCP pe-STARR-Seq --------------------------------------------------------------------#
  # Sanger Seq
  file.edit("git_peSTARRSeq/subscripts/sanger_peSTARRSeq_SCR1_luc_validations.R")
  # Process validations data
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
  # Compare o/e (residuals) for STARR-Seq vs luciferase
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_vs_STARRSeq_residuals.R")
  # plot VALIDATIONS
  file.edit("git_peSTARRSeq/subscripts/PCC_STARRSeq_luciferase_validations.R")
  file.edit("git_peSTARRSeq/subscripts/barplot_luciferase_validations.R")
  file.edit("git_peSTARRSeq/subscripts/connected_stripchart_luciferase_validations.R")
}