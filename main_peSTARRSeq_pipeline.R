setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

#----------------------------------------------------------------------------------------------------#
# peSTARR-Seq pipeline
#----------------------------------------------------------------------------------------------------#
if(F)
{
  # Used by the pipeline (indexes and functions)
  file.edit("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/create_twist8_subread_index.R")
  file.edit("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/create_twist12_subread_index.R")
  #------ Pipeline 3.0 ------#
  file.edit("git_peSTARRSeq/functions/pipeline_3.0.R")
  file.edit("git_peSTARRSeq/subscripts/run_pipeline_3.0_parallel.R") # Parallel!
  file.edit("git_peSTARRSeq/subscripts/rename_files_pipeline_3.0.R")
  file.edit("git_peSTARRSeq/subscripts/multi_screen_models.R")
  #----------- QC -----------#
  source("git_peSTARRSeq/subscripts/aggregate_alignment_statistics.R")
  source("git_peSTARRSeq/subscripts/sequencing_saturation.R")
  source("git_peSTARRSeq/subscripts/barplot_read_per_theoretical_pair.R")
  source("git_peSTARRSeq/subscripts/PCC.R")
}
