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
  # Hamlet luciferase assays
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/HAM_luciferase_assays.R")
  
  #-------------------------#
  # rev/intron STARRSeq design
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/enh_promoter_distance.R")
  file.edit("git_peSTARRSeq/subscripts/intron_selection.R")
  file.edit("git_peSTARRSeq/subscripts/design_introns_primers.R")
  
  file.edit("git_peSTARRSeq/subscripts/restriction_enzyme_selection.R")
  
  
  #-------------------------#
  # peSTARR-Seq pipeline
  #-------------------------#

}
