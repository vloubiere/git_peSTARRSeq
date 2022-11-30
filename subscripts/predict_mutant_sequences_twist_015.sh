#!/bin/bash
cd /groups/stark/vloubiere/projects/pe_STARRSeq/db/STARR_predictions/twist_lib_015/
FASTA=/groups/stark/vloubiere/projects/pe_STARRSeq/db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_with_WT.fasta
bsub=/groups/stark/software-all/shell/bsub_gridengine

### Split multiple FASTA
GENOME=dm3
$bsub -o log_fasta "module load build-env/2020; module load kent_tools/20190507-linux.x86_64; faSplit sequence $FASTA 500 FASTA_seqs/Split_"

#---------------------------------#
# Wait previous chunk is completed before running this one
#---------------------------------#
### predict DeepSTARR scores
pred_script=/groups/stark/almeida/Projects/High_resolution_STARRseq/results/20210122_deep_learning_prediction_bins_100bp_stride/Predict_dev_and_hk_enh_activity_CNN_model_from_fasta.py
Model_name=/groups/stark/almeida/Projects/High_resolution_STARRseq/results/20210122_deep_learning_prediction_bins_100bp_stride/Final_CNN_models/Model_continuous_CNN_pool_120_7_60_3_60_3_pooling_2denses_dropout0.3_valloss2.73

for file in `ls FASTA_seqs/*fa`; do
  $bsub -o log_predict -n predict -T 20:00 "$pred_script -s $file -m $Model_name"
done

### Combine results
head -n 1 FASTA_seqs/Split_000.fa_predictions_Model_continuous_CNN_pool_120_7_60_3_60_3_pooling_2denses_dropout0.3_valloss2.73.txt > twist_0015_mutated_sequences_with_WT.fasta_DeepSTARR_predictions.txt
tail -n +2 -q FASTA_seqs/*.txt >> twist_0015_mutated_sequences_with_WT.fasta_DeepSTARR_predictions.txt