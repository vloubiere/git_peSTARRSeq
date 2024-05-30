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
# Luciferase processing ----
file.edit("git_peSTARRSeq/subscripts/validations_luciferase_data_processing.R")
file.edit("git_peSTARRSeq/subscripts/homotypic_validations_data_processing.R")
# Linear models ----
file.edit("git_peSTARRSeq/subscripts/linear_model_large_WT_lib.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_longer_spacer_act_pairs.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_ECD_full_dataset.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_OSC_full_dataset.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_focused_RpS12_lib.R")
file.edit("git_peSTARRSeq/subscripts/linear_model_focused_DSCP_lib.R")
# Motif analyses ----
file.edit("git_peSTARRSeq/subscripts/Motif_counts_large_WT_oligo_pool.R")
file.edit("git_peSTARRSeq/subscripts/Motif_combinations_compute_mean_act_and_residuals.R")
# LASSO models ----
file.edit("git_peSTARRSeq/subscripts/lasso_model_DSCP_large_WT.R")
file.edit("git_peSTARRSeq/subscripts/lasso_model_ECD.R")
file.edit("git_peSTARRSeq/subscripts/lasso_model_OSC.R")
# RF models ----
file.edit("git_peSTARRSeq/subscripts/randomForest_model_DSCP_large_WT.R")
file.edit("git_peSTARRSeq/subscripts/randomForest_model_ECD_WT.R")
file.edit("git_peSTARRSeq/subscripts/randomForest_model_OSC_WT.R")
# Define developmental/Hk genes for IDR content ----
file.edit("git_peSTARRSeq/subscripts/define_dev_hk_genes.R")

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
file.edit("git_peSTARRSeq/subscripts/activity_dist_enhancer_pairs_activity_matched_controls.R")
pdftools::pdf_combine(c("pdf/draft/histogram_residuals.pdf",
                        "pdf/draft/motif_impact_act_vs_res.pdf",
                        "pdf/draft/mutant_library_per_cdition.pdf"),
                      "pdf/figures/Figure_3.pdf")
file.show("pdf/figures/Figure_3.pdf")

## Figure 4 ----
file.edit("git_peSTARRSeq/subscripts/hk_additive_vs_multiplicative_model.R")
file.edit("git_peSTARRSeq/subscripts/hk_accuracy_predictions_pie_charts.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_enh_vs_CP.R")
file.edit("git_peSTARRSeq/subscripts/review_examples_additive_multiplicative.R")
file.edit("git_peSTARRSeq/subscripts/hk_vs_dev_COFs_TFs_IDR_content.R")
pdftools::pdf_combine(c("pdf/draft/Modelling_obs_vs_expected_focused_RpS12_lib.pdf",
                        "pdf/draft/Hk_accuracy_predictions_pie_charts.pdf",
                        "pdf/draft/boxplot_hkCP_dev_vs_hk_pairs.pdf",
                        "pdf/draft/review_examples_add_mult_hkCP_dCP.pdf",
                        "pdf/draft/hk_dev_TFs_COFs_IDR_content.pdf"),
                      "pdf/figures/Figure_4.pdf")
file.show("pdf/figures/Figure_4.pdf")

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
file.edit("git_peSTARRSeq/subscripts/enh_enh_dist_genome.R")
file.edit("git_peSTARRSeq/subscripts/model_longer_spacer.R")
file.edit("git_peSTARRSeq/subscripts/ecd_OSC_induction.R")
file.edit("git_peSTARRSeq/subscripts/luciferase_homotypic_pairs.R")
pdftools::pdf_combine(c("pdf/draft/Compare_rsq_models.pdf",
                        "pdf/draft/Accuracy_predictions_pie_charts_with_lm.pdf",
                        "pdf/draft/Compare_mult_linear_predictions_vllib002.pdf",# Quantification of additive residuals plotted in Fig 2
                        "pdf/draft/promoter_saturation.pdf",
                        "pdf/draft/ecd_OSC_induction.pdf",
                        "pdf/draft/super_additivity_OSC_ECD.pdf"),
                      "pdf/figures/Supplementary_figure_2.pdf")
file.show("pdf/figures/Supplementary_figure_2.pdf")

# Supplementary Figure 3 ----
file.edit("git_peSTARRSeq/subscripts/Motif_combinations_homotypic_vs_hetero.R")
file.edit("git_peSTARRSeq/subscripts/evaluate_lasso_model_activity.R")
file.edit("git_peSTARRSeq/subscripts/LASSO_coeffs_per_cell_type.R")
file.edit("git_peSTARRSeq/subscripts/evaluate_lasso_model_residuals.R")
file.edit("git_peSTARRSeq/subscripts/hk_dev_pairs_vs_CP.R")
file.edit("git_peSTARRSeq/subscripts/motifs_activity_single_enhancers.R")
pdftools::pdf_combine(c("pdf/draft/motif_homotypic_vs_hetero_act_vs_res.pdf",
                        "pdf/draft/lasso_large_WT_activity_prediction.pdf",
                        "pdf/draft/LASSO_coeff_per_cell_type.pdf",
                        "pdf/draft/lasso_residuals_per_cdition.pdf",
                        "pdf/draft/mutant_library_per_cdition.pdf"),# Impact of Dref motif plotted in fig 3
                      "pdf/figures/Supplementary_figure_3.pdf")
file.show("pdf/figures/Supplementary_figure_3.pdf")

# Supplementary tables
file.edit("git_peSTARRSeq/subscripts/supplementary_tables.R")

# GEO submission ----
file.edit("git_peSTARRSeq/subscripts/GEO_submission.R")
file.edit("git_peSTARRSeq/subscripts/GEO_submission_review.R")

# Figures to reviewers ----
file.edit("git_peSTARRSeq/subscripts/Motif_combinations_homotypic_vs_hetero_test.R")
file.edit("git_peSTARRSeq/subscripts/review_orientation_luciferase.R") # Luciferase enh/ctl ctl/enh
file.edit("git_peSTARRSeq/subscripts/review_PROSeq_hk_dev.R") # Compare housekeeping developmental PROSeq
file.edit("git_peSTARRSeq/subscripts/mutant_library_wt_residuals.R")
pdftools::pdf_combine(c("pdf/draft/review_luciferase_vs_STARRSeq_differences.pdf",
                        "pdf/draft/review_PROSeq_hk_vs_dev_genes.pdf",
                        "pdf/draft/reviewers_long_spacer.pdf"),
                      "pdf/figures/Figures_to_reviewers.pdf")
file.show("pdf/figures/Figures_to_reviewers.pdf")
file.edit("git_peSTARRSeq/subscripts/analysis_mutant_lib.R")

# Review ----
file.edit("git_peSTARRSeq/subscripts/shared_enhancers_hkCP.R")
file.edit("git_peSTARRSeq/subscripts/longer_spacer_individual_act.R")
# Luciferase super-additivity
file.edit("git_peSTARRSeq/subscripts/review_luciferase_super_additivity.R")
# Trl twist motifs residuals
file.edit("git_peSTARRSeq/subscripts/review_trl_twist_residuals.R")
file.edit("git_peSTARRSeq/subscripts/review_trl_twist_residuals_boxplot.R")
# Dev Hk heterotypic pairs
file.edit("git_peSTARRSeq/subscripts/review_dev_hk_heterotypic pairs.R")
# ECD/OSC screens
file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_ECD.R")
file.edit("git_peSTARRSeq/subscripts/Compare_additive_multiplicative_OSC.R")
# Homotypic enhancer pairs
file.edit("git_peSTARRSeq/subscripts/ham_luciferase_data_processing.R") # Luciferase homotypic pairs
file.edit("git_peSTARRSeq/subscripts/review_ham_enhancer_homotypic_pairs_luciferase.R")
# Downsample WT lib ----
file.edit("git_peSTARRSeq/subscripts/downsample_large_WT_lib.R") # Downsample WT library -> affects result?
