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
  file.edit("git_peSTARRSeq/subscripts/activity_aggregate.R")
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