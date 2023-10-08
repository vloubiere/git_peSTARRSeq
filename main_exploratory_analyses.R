setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# QC and sanity check ----
file.edit("git_peSTARRSeq/subscripts/Benchmarking_UMI_collapsing.R") # Compare to BA approach
file.edit("git_peSTARRSeq/subscripts/Compare_DESeq2_to_raw_log2_ratio.R") # Sanity check

# Objects ----
file.edit("git_peSTARRSeq/subscripts/PCA_vllib002_residuals.R")
file.edit("git_peSTARRSeq/subscripts/PCA_vllib002_actPairs_residuals.R")

# Alignment statistics ----
file.edit("git_peSTARRSeq/subscripts/alignment_statistics.R")

# Candidates classification ----
file.edit("git_peSTARRSeq/subscripts/Candidates_classification_vllib002.R")

# Modelling vllib002  ----
file.edit("git_peSTARRSeq/subscripts/scheme_add_mult.R")
file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_vllib002_actPairs.R")

# Residuals modelling and motifs impact ----
## Residuals heatmap ordered by activity ----
file.edit("git_peSTARRSeq/subscripts/heatmap_res_ordered_ind_act_actBreaks_vllib002.R")
file.edit("git_peSTARRSeq/subscripts/heatmap_res_ordered_ind_act_vllib002.R")
file.edit("git_peSTARRSeq/subscripts/heatmap_res_ordered_ind_act_vllib002_classes.R")
file.edit("git_peSTARRSeq/subscripts/heatmap_res_ordered_ind_act_inact_pairs_vllib002.R")
## Motifs impact vs activity matched controls ----
file.edit("git_peSTARRSeq/subscripts/activity_predictive_motifs_activity_matched_controls.R")
## Residuals heatmap ordered by residuals ----
file.edit("git_peSTARRSeq/subscripts/heatmap_residuals_ordered_pca_vllib002.R")
file.edit("git_peSTARRSeq/subscripts/heatmap_residuals_ordered_pca_actPairs_vllib002.R")

# Enhancers from same genes vs diff genes ----
file.edit("git_peSTARRSeq/subscripts/activity_shn_enhancer_pairs_activity_matched_controls.R")

# Hk vs dev ----
file.edit("git_peSTARRSeq/subscripts/hk_motif_impact_act_vs_res.R")

# longer spacer ----
file.edit("git_peSTARRSeq/subscripts/enh_enh_dist_genome.R")
file.edit("git_peSTARRSeq/subscripts/long_spacer.R")

# Alterative promoters ----
file.edit("git_peSTARRSeq/subscripts/alternative_promoters.R")

# DHS library ----
file.edit("git_peSTARRSeq/subscripts/linear_model_vllib030.R") # train linear model on the whole library
file.edit("git_peSTARRSeq/subscripts/DHS_library.R")

# Van Steensel data ----
file.edit("git_peSTARRSeq/subscripts/analyse_vanSteensel_data.R")

# Wrappers to clean folder ----
file.edit("git_peSTARRSeq/subscripts/clean_scripts_pdf.R")

# Old/not used ----
file.edit("git_peSTARRSeq/subscripts/motif_enrichment_residuals.R")
# Saturation ----
file.edit("git_peSTARRSeq/subscripts/umi_saturation.R") # No saturation of the UMI
file.edit("git_peSTARRSeq/subscripts/dev_enhancer_stength_vs_cooperativity.R") # Strong enhancers dont synergize and dont show specific motifs
## Analysis vllib002 residuals Trl/Twist ----
file.edit("git_peSTARRSeq/subscripts/lm_residuals_examples_dev_pairs_vllib002.R") # Linear model examples
file.edit("git_peSTARRSeq/subscripts/density_residuals_motif_enrichment_dev_pairs_vllib002.R") # Show activity/residuals relationship
file.edit("git_peSTARRSeq/subscripts/heatmaps_motifs_residuals_dev_pairs_vllib002.R") # Heatmap lm coeffs to predict residuals
file.edit("git_peSTARRSeq/subscripts/mutant_library_sequencing_depth.R")

