setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Luciferase validations ----
if(F)
{
  # HAM enhancer DSCP (Lorena) ----
  ## Sanger seq ----
  file.edit("git_peSTARRSeq/subscripts/sanger_HAM_SUP1_luciferase_constructs.R")
  ## Process validations data ----
  file.edit("git_peSTARRSeq/subscripts/ham_luciferase_data_processing.R")
  ## Luciferase ----
  file.edit("git_peSTARRSeq/subscripts/HAM_luciferase_assays.R") # Luciferase
  
  # dev enhancers DSCP pe-STARR-Seq large validation ----
  ## Sanger Seq ----
  file.edit("git_peSTARRSeq/subscripts/sanger_peSTARRSeq_SCR1_luc_validations.R")
  ## Process validations data ----
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
  ## Compare o/e (residuals) for STARR-Seq vs luciferase ----
  file.edit("git_peSTARRSeq/subscripts/validations_luciferase_vs_STARRSeq_residuals.R")
  ## plot VALIDATIONS ----
  file.edit("git_peSTARRSeq/subscripts/PCC_STARRSeq_luciferase_validations.R")
  file.edit("git_peSTARRSeq/subscripts/barplot_luciferase_validations.R")
  file.edit("git_peSTARRSeq/subscripts/connected_stripchart_luciferase_validations.R")
  
  # Paper review homotypic enhancer pairs ----
  file.edit("git_peSTARRSeq/subscripts/review_homotypic_enhancer_pairs_PCR_design.R")
  file.edit("git_peSTARRSeq/subscripts/review_homotypic_enhancer_pairs_sequences.R")
  
  file.edit("git_peSTARRSeq/subscripts/homotypic_validations_data_processing.R")
  
}