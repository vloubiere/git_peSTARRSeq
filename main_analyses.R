setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

if(F)
{
  #-------------------------#
  # Final tables
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/make_300bp_uniq_enhancer_features_object.R")
  file.edit("git_peSTARRSeq/functions/generate_final_data_table.R")
  
  #-------------------------#
  # Intron impact
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/PCC_input_w_wo_intron_INPUT.R")
  file.edit("git_peSTARRSeq/subscripts/PCC_input_w_wo_intron_SCREEN.R")
  
  #-------------------------#
  # Observed vs additive
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/Figure_1AB.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_1CD.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_1E.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_1FG.R")
  
  file.edit("git_peSTARRSeq/subscripts/Figure_2A.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_2B.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_2C.R")  
  file.edit("git_peSTARRSeq/subscripts/Figure_2D.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_2E.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_2F.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_2G.R")
  
  file.edit("git_peSTARRSeq/subscripts/Figure_3AB.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_3C.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_3D.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_3E.R")
  
  file.edit("git_peSTARRSeq/subscripts/Figure_4A.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_4B.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_4C.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_4D.R")
  
  #-------------------------#
  # Supplementary figures
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/enh_enh_distance.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_inact_enh_pairs_heatmap.R")
  
  #-------------------------#
  # Exploration
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_differences.R")
  file.edit("git_peSTARRSeq/subscripts/gene_activity_prediction.R")
  file.edit("git_peSTARRSeq/subscripts/Sup_Figure_intron_spacers.R")
  file.edit("git_peSTARRSeq/subscripts/Sup_sequencing_depth_vs_robustness.R")
  file.edit("git_peSTARRSeq/subscripts/spearman_residuals_marginals.R")
  
  #-------------------------#
  # Tests
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/activity_aggregate.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_DSCP_large_screen.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_RpS12_screen.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_compare_hk_dev.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_RpS12_inactive_pairs.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_RpS12_DHS_peaks.R")
  file.edit("git_peSTARRSeq/subscripts/Figure_DSCP_screen.R")
  file.edit("git_peSTARRSeq/subscripts/lasso_modelling_dev_hk.R")
  file.edit("git_peSTARRSeq/subscripts/lm_modelling_dev_hk.R")
  file.edit("git_peSTARRSeq/subscripts/randomForest_dev_hk.R")
  file.edit("git_peSTARRSeq/subscripts/feature_selection.R")
  file.edit("git_peSTARRSeq/subscripts/lm_modelling_dev_hk.R")
  file.edit("git_peSTARRSeq/subscripts/smoothScatters_additivity_dev_vs_hk.R")
  file.edit("test/tests.R")
  file.edit("test/tests.R")
  file.edit("test/tests2.R")
}