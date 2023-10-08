setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# pSTARR-Seq PROCESSING ----
## Functions ----
file.edit("git_peSTARRSeq/functions/pSTARRSeq_pipeline.R")
file.edit("git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R")
## Alignment indexes ----
file.edit("git_peSTARRSeq/subscripts/create_twist8_subread_index.R")
file.edit("git_peSTARRSeq/subscripts/create_twist12_subread_index.R")
file.edit("git_peSTARRSeq/subscripts/create_twist15_subread_index.R")
## Run pipeline ----
file.edit("git_peSTARRSeq/subscripts/extract_reads_VBC_bam.R")
file.edit("git_peSTARRSeq/subscripts/run_pipeline_parallel.R")
file.edit("git_peSTARRSeq/subscripts/vllib029_vllib030_pseudo_reps.R")
file.edit("git_peSTARRSeq/subscripts/compute_log2FoldChange.R")

# Figures ----
## Figure 1 ----
file.edit("git_peSTARRSeq/subscripts/Correlation_left_rigth_activities.R")
file.edit("git_peSTARRSeq/subscripts/Compare_individual_vs_enh_pairs.R")
file.edit("git_peSTARRSeq/subscripts/heatmap_act_ordered_ind_act_vllib002.R")

## Figure 2 ----
file.edit("git_peSTARRSeq/subscripts/linear_model_vllib002.R")
file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_vllib002.R")

## Figure 3 ----
file.edit("git_peSTARRSeq/subscripts/Motif_counts_large_WT_oligo_pool.R")
file.edit("git_peSTARRSeq/subscripts/motif_impact_act_vs_res.R")
file.edit("git_peSTARRSeq/subscripts/mutant_library.R")
file.edit("git_peSTARRSeq/subscripts/activity_dist_enhancer_pairs_activity_matched_controls.R")

## Figure 4 ----
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_additive_vs_linear_model.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_enh_vs_CP.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_IDR_fraction.R")
file.edit("git_peSTARRSeq/subscripts/define_dev_hk_genes.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_clustering.R")

# Supplementary Figure 1 ----
file.edit("git_peSTARRSeq/functions/characterization_libraries.R")
file.edit("git_peSTARRSeq/subscripts/PCC_peSTARRSeq_replicates.R")
file.edit("git_peSTARRSeq/subscripts/Luciferase_validations.R")

# Supplementary Figure 2 ----
# Rsqaured barplot (comparing the different models) and accuracy pie chart are included in Figure 2
file.edit("git_peSTARRSeq/subscripts/Compare_multiplicative_linear_predicted.R")
file.edit("git_peSTARRSeq/subscripts/promoter_saturation.R")

# Supplementary Figure 3 ----
file.edit("git_peSTARRSeq/subscripts/lasso_model_vllib002.R")
file.edit("git_peSTARRSeq/subscripts/evaluate_lasso_model.R")

# Supplementary tables
file.edit("git_peSTARRSeq/subscripts/supplementary_tables.R")

# GEO submission ----
file.edit("git_peSTARRSeq/subscripts/GEO_submission.R")
