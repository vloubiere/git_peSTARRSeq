setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

if(F)
{
  #-------------------------#
  # HAM Luciferase assays
  #-------------------------#
  source("git_peSTARRSeq/subscripts/HAM_luciferase_assays.R")
  
  #-------------------------#
  # STARR-Seq VALIDATIONS
  #-------------------------#
  # Process validations data
  source("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
  # Compare o/e (residuals) for STARR-Seq vs luciferase
  source("git_peSTARRSeq/subscripts/validations_luciferase_vs_STARRSeq_residuals.R")
  # plot VALIDATIONS
  source("git_peSTARRSeq/subscripts/PCC_STARRSeq_luciferase_validations.R")
  source("git_peSTARRSeq/subscripts/barplot_luciferase_validations.R")
  source("git_peSTARRSeq/subscripts/connected_stripchart_luciferase_validations.R")
  
  #-------------------------#
  # Sanger sequencing
  #-------------------------#
  # Hamlet
  file.edit("git_peSTARRSeq/subscripts/sanger_HAM_SUP1_luciferase_constructs.R")
  # Validations
  file.edit("git_peSTARRSeq/subscripts/sanger_peSTARRSeq_SCR1_luc_validations.R")
}
