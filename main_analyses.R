setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

if(F)
{
  #-------------------------#
  # Final tables
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/make_300bp_uniq_enhancer_features_object.R")
  file.edit("git_peSTARRSeq/functions/generate_final_data_table.R")
  
  #-------------------------#
  # Aggregate activity
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/activity_aggregate.R")
  file.edit("git_peSTARRSeq/subscripts/smoothsScatter_obs_expAdd.R")
  file.edit("git_peSTARRSeq/subscripts/boxplot_FC_silencer_insulator_LR.R")
  
  #-------------------------#
  # Fine grained heatmaps
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/heatmap_activity.R")
  file.edit("git_peSTARRSeq/subscripts/detailed_heatmap_additivity.R")
  
  #-------------------------#
  # peSTARR-Seq modeling
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/full_data_linear_models_activity_residuals.R")
  file.edit("git_peSTARRSeq/subscripts/collapsed_data_linear_models_activity_residuals.R")
  file.edit("git_peSTARRSeq/subscripts/smoothScatterplot_motifs_contribution.R")
  file.edit("git_peSTARRSeq/subscripts/barplot_motifs_contribution.R")
  
  #-------------------------#
  # Modelling gene endogenous activity
  #-------------------------#
  file.edit("git_peSTARRSeq/subscripts/Modelling_endogenous_activity.R")
  
  #-------------------------#
  # hk/dev specificity in paired STARR-Seq
  #-------------------------#
  
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
}
