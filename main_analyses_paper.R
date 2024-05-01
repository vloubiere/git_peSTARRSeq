setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# pSTARR-Seq PROCESSING ----
## Functions ----
file.edit("git_peSTARRSeq/functions/pSTARRSeq_pipeline.R")
file.edit("git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R")
file.edit("git_peSTARRSeq/functions/train_LASSO_model.R")
file.edit("git_peSTARRSeq/functions/train_rf_model.R")
## Alignment indexes ----
file.edit("git_peSTARRSeq/subscripts/create_twist8_subread_index.R")
file.edit("git_peSTARRSeq/subscripts/create_twist12_subread_index.R")
file.edit("git_peSTARRSeq/subscripts/create_twist15_subread_index.R")
## Run pipeline ----
file.edit("git_peSTARRSeq/subscripts/extract_reads_VBC_bam.R")
file.edit("git_peSTARRSeq/subscripts/run_pipeline_parallel.R")
file.edit("git_peSTARRSeq/subscripts/compute_log2FoldChange.R")
## Downstream analyses and tables ----
# Linear models ----
file.edit("git_peSTARRSeq/subscripts/linear_model_large_WT_lib.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_ECD.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_OSC.R")
# Motif analyses ----
file.edit("git_peSTARRSeq/subscripts/Motif_counts_large_WT_oligo_pool.R")
file.edit("git_peSTARRSeq/subscripts/Motif_combinations_compute_mean_act_and_residuals.R")
# LASSO models ----
file.edit("git_peSTARRSeq/subscripts/lasso_model_DSCP_large_WT.R")
file.edit("git_peSTARRSeq/subscripts/lasso_model_ECD_WT.R")
file.edit("git_peSTARRSeq/subscripts/lasso_model_OSC_WT.R")
# RF models ----
file.edit("git_peSTARRSeq/subscripts/randomForest_model_DSCP_large_WT.R")
file.edit("git_peSTARRSeq/subscripts/randomForest_model_ECD_WT.R")
file.edit("git_peSTARRSeq/subscripts/randomForest_model_OSC_WT.R")

# Figures ----
## Figure 1 ----
file.edit("git_peSTARRSeq/subscripts/Correlation_ind_left_rigth_activities.R")
file.edit("git_peSTARRSeq/subscripts/Compare_individual_vs_enh_pairs.R")
file.edit("git_peSTARRSeq/subscripts/heatmap_act_ordered_ind_act_vllib002.R")
file.edit("git_peSTARRSeq/subscripts/review_pair_orientation_vs_act.R")
pdftools::pdf_combine(c("pdf/draft/Correlation_left_rigth_activities.pdf",
                        "pdf/draft/Compare_individual_vs_enh_pairs.pdf",
                        "pdf/draft/heatmap_ordered_ind_act.pdf",
                        "pdf/draft/review_orientation_activity.pdf"),
                      "pdf/figures/Figure_1.pdf")
file.show("pdf/figures/Figure_1.pdf")

## Figure 2 ----
file.edit("git_peSTARRSeq/subscripts/models_obs_vs_exp_large_WT_lib.R")
file.edit("git_peSTARRSeq/subscripts/accuracy_predictions_pie_charts.R")
file.edit("git_peSTARRSeq/subscripts/models_ECD_OSC.R")
pdftools::pdf_combine(c("pdf/draft/Modelling_obs_vs_expected_large_WT_lib.pdf",
                        "pdf/draft/Accuracy_predictions_pie_charts.pdf",
                        "pdf/draft/super_additivity_OSC_ECD.pdf"),
                      "pdf/figures/Figure_2.pdf")
file.show("pdf/figures/Figure_2.pdf")

## Figure 3 ----
file.edit("git_peSTARRSeq/subscripts/histogram_residuals_large_WT_lib.R")
file.edit("git_peSTARRSeq/subscripts/Motif_combinations_compare_mean_act_vs_res.R")
file.edit("git_peSTARRSeq/subscripts/mutant_library.R")
pdftools::pdf_combine(c("pdf/draft/histogram_residuals.pdf",
                        "pdf/draft/motif_impact_act_vs_res.pdf",
                        "pdf/draft/mutant_library_per_cdition.pdf"),
                      "pdf/figures/Figure_3.pdf")
file.show("pdf/figures/Figure_3.pdf")
file.edit("git_peSTARRSeq/subscripts/activity_dist_enhancer_pairs_activity_matched_controls.R") # Removed following reviewer's suggestion

## Figure 4 ----
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_additive_vs_linear_model.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_enh_vs_CP.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_IDR_fraction.R")
file.edit("git_peSTARRSeq/subscripts/define_dev_hk_genes.R")
file.edit("git_peSTARRSeq/subscripts/review_examples_additive_multiplicative.R") # Additive vs multiplicative examples
file.edit("git_peSTARRSeq/subscripts/review_hk_vs_dev_TF_length.R") # Hk/Dev TFs length

# Supplementary Figure 1 ----
file.edit("git_peSTARRSeq/subscripts/PCC_peSTARRSeq_replicates.R")
file.edit("git_peSTARRSeq/subscripts/Luciferase_validations.R")
pdftools::pdf_combine(c("pdf/draft/Luciferase_validations.pdf",
                        "pdf/draft/PCC_peSTARRSeq_replicates.pdf"),
                      "pdf/figures/Supplementary_figure_1.pdf")
file.show("pdf/figures/Supplementary_figure_1.pdf")

# Supplementary Figure 2 ----
file.edit("git_peSTARRSeq/subscripts/barplot_rsq_compare_models.R")
file.edit("git_peSTARRSeq/subscripts/accuracy_pie_charts_with_lm.R")
file.edit("git_peSTARRSeq/subscripts/Compare_multiplicative_linear_predicted.R")
file.edit("git_peSTARRSeq/subscripts/promoter_saturation.R")
# Quantification of additive residuals plotted in Fig 2
file.edit("git_peSTARRSeq/subscripts/ecd_OSC_induction.R")
pdftools::pdf_combine(c("pdf/draft/Compare_rsq_models.pdf",
                        "pdf/draft/Accuracy_predictions_pie_charts_with_lm.pdf",
                        "pdf/draft/Compare_mult_linear_predictions_vllib002.pdf",
                        "pdf/draft/promoter_saturation.pdf",
                        "pdf/draft/ecd_OSC_induction.pdf",
                        "pdf/draft/super_additivity_OSC_ECD.pdf"),
                      "pdf/figures/Supplementary_figure_2.pdf")
file.show("pdf/figures/Supplementary_figure_2.pdf")

# Supplementary Figure 3 ----
file.edit("git_peSTARRSeq/subscripts/Motif_combinations_homotypic_vs_hetero.R")
file.edit("git_peSTARRSeq/subscripts/evaluate_lasso_model.R")
file.edit("git_peSTARRSeq/subscripts/LASSO_coeffs_per_cell_type.R")
pdftools::pdf_combine(c("pdf/draft/motif_homotypic_vs_hetero_act_vs_res.pdf",
                        "pdf/draft/residuals_prediction_lasso.pdf"),
                      "pdf/figures/Supplementary_figure_3.pdf")
file.show("pdf/figures/Supplementary_figure_3.pdf")


# Supplementary tables
file.edit("git_peSTARRSeq/subscripts/supplementary_tables.R")

# GEO submission ----
file.edit("git_peSTARRSeq/subscripts/GEO_submission.R")

# Review ----
# Luciferase super-additivity
file.edit("git_peSTARRSeq/subscripts/review_luciferase_super_additivity.R")
# Ecdysone screen
file.edit("git_peSTARRSeq/subscripts/review_ECD_OSC_enhancer_pairs.R")
file.edit("git_peSTARRSeq/subscripts/review_ECD_OSC_coop.R")
# Trl twist motifs residuals
file.edit("git_peSTARRSeq/subscripts/review_trl_twist_residuals.R")
file.edit("git_peSTARRSeq/subscripts/review_trl_twist_residuals_boxplot.R")
# Dev Hk heterotypic pairs
file.edit("git_peSTARRSeq/subscripts/review_dev_hk_heterotypic pairs.R")
# ECD screen
file.edit("git_peSTARRSeq/subscripts/linear_model_ECD.R")
file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_ECD.R")
# OSC screen
file.edit("git_peSTARRSeq/subscripts/linear_model_OSC.R")
file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_OSC.R")

# Homotypic enhancer pairs
file.edit("git_peSTARRSeq/subscripts/ham_luciferase_data_processing.R") # Luciferase homotypic pairs
file.edit("git_peSTARRSeq/subscripts/review_homotypic_luciferase.R")
file.edit("git_peSTARRSeq/subscripts/review_homotypic_PCR_design.R")

# Figures to reviewers ----
file.edit("git_peSTARRSeq/subscripts/review_orientation_luciferase.R") # Luciferase enh/ctl ctl/enh
file.edit("git_peSTARRSeq/subscripts/review_PROSeq_hk_dev.R") # Compare housekeeping developmental PROSeq
file.edit("git_peSTARRSeq/subscripts/mutant_library_wt_residuals.R")
pdftools::pdf_combine(c("pdf/draft/review_luciferase_vs_STARRSeq_differences.pdf",
                        "pdf/draft/review_PROSeq_hk_vs_dev_genes.pdf"),
                      "pdf/figures/Figures_to_reviewers.pdf")
file.show("pdf/figures/Figures_to_reviewers.pdf")

# Tests/unused ----
file.edit("git_peSTARRSeq/subscripts/review_synthetic_data.R") # Create add/mult synthetic datasets
file.edit("git_peSTARRSeq/subscripts/downsample_large_WT_lib.R") # Downsample WT library -> affects result?
file.edit("git_peSTARRSeq/subscripts/tests_mutant_library.R")