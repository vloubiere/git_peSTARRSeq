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
  file.edit("git_peSTARRSeq/subscripts/SUHW_peak_selection.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_peak_calling.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_enrichment_quantif.R")
  file.edit("git_peSTARRSeq/subscripts/ATAC_STARR_DHS_candidates_selection.R")# Somehow, peaks are different from original run. See old folders
  file.edit("git_peSTARRSeq/subscripts/CPs_selection_twist.R")
  file.edit("git_peSTARRSeq/subscripts/design_twist_012.R") # lib can be found at "Rdata/vl_library_twist12_210610.rds"
  
  #-------------------------#
  # STARRSeq design
  #-------------------------#
  # revSTARRSeq
  file.edit("git_peSTARRSeq/subscripts/restriction_enzyme_selection.R") # NotI cuts the peSTARRSeq vector but is not detected!?
  file.edit("git_peSTARRSeq/subscripts/spacer_selection_revSTARRSeq.R")
  file.edit("git_peSTARRSeq/subscripts/primer_design_revSTARRSeq_spacers.R")
  file.edit("git_peSTARRSeq/subscripts/enh_promoter_distance.R")
  file.edit("git_peSTARRSeq/subscripts/qPCR_non_self_lig_vllib004.R")
  # Introns without enhancer
  file.edit("git_peSTARRSeq/subscripts/intron_selection.R")
  file.edit("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/design_introns_primers.R")
  # Introns with enhancer
  file.edit("git_peSTARRSeq/subscripts/intron_enhancer_selection.R")
  file.edit("git_peSTARRSeq/subscripts/design_introns_enhancer_primers.R")
  # Change promoter
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_to_replace_DSCP.R") # Forgot to check that CPs do not contain AgeI/SalI sites (see next line)
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_no_enzymatic_restriciton_to_replace_DSCP.R") # I realized the design was not optimal cause did not take into account hk induction!
  file.edit("git_peSTARRSeq/subscripts/select_active_CP_to_replace_DSCP_final.R") # I realized the design was not optimal cause did not take into account hk induction!
  
  #-------------------------#
  # Sanger sequencing
  #-------------------------#
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
