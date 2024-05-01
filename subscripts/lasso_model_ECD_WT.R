setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)

# Predict values for the full dataset (not only ecd pairs)
if(!file.exists("db/linear_models/FC_DSCP_ECD_WT_lm_predictions_full_data_LASSO.rds"))
{
  # Import data and remove control pairs ----
  dat <- readRDS("db/FC_tables/DSCP_ECD_WT_FC_DESeq2.rds")
  dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
  
  # Predicted values ----
  dat[, `Additive model`:= log2(2^indL+2^indR-1)]
  dat[, `Multiplicative model`:= indL+indR]
  model <- lm(log2FoldChange~indL*indR, dat)
  dat[, `Linear model`:= predict(model)]
  dat[, residuals:= log2FoldChange-`Linear model`]
  
  # Define train and test sets ----
  set.seed(1)
  dat[dat[, .(set= sample(3)), L], setL:= i.set, on= "L"]
  set.seed(1)
  dat[dat[, .(set= sample(3)), R], setR:= i.set, on= "R"]
  dat[, set:= .GRP, .(setL, setR)]
  dat$setL <- dat$setR <- NULL
  
  saveRDS(dat,
          "db/linear_models/FC_DSCP_ECD_WT_lm_predictions_full_data_LASSO.rds")
}else
  dat <- readRDS("db/linear_models/FC_DSCP_ECD_WT_lm_predictions_full_data_LASSO.rds")

# Full dataset (no pairs exclude, all pairs predicted) ----
fullData <- data.table(set= 0,
                       excludeL= "None", # Use all left enhancers
                       excludeR= "None") # Use all right enhancers

# Define 10 folds to exclude from training data and CV ----
set <- dat[, {
  .(excludeL= paste0(unique(L), collapse = ","),
    excludeR= paste0(unique(R), collapse = ","))
}, set]
set <- rbind(set, fullData)

# Predict act and residuals ----
set[, output:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/lasso_models/",
                      ifelse(set==0, "full_dataset", paste0("fold", set)), "_ECD_WT_lasso.rds")]

# Command ----
set[, cmd:= {
  if(!file.exists(output))
  {
    # [required] 1/ Comma-separated list of left enhancers to exclude from training\n
    # [required] 2/ Comma-separated list of right enhancers to exclude from training\n
    # [required] 3/ Motif counts matrix with first column containing enhancer IDs \n
    # [required] 4/ .rds input data file containing the response variables \n"
    # [required] 5/ .rds output file \n
    paste("/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
          "/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/train_LASSO_model.R",
          excludeL,
          excludeR,
          "/groups/stark/vloubiere/projects/pe_STARRSeq/db/motif_counts/twist008_motif_counts_selected.rds",
          "/groups/stark/vloubiere/projects/pe_STARRSeq/db/linear_models/FC_DSCP_ECD_WT_lm_predictions_full_data_LASSO.rds",
          output)
  }
}, output]

# Submits ----
run <- set[!is.na(cmd)]
if(nrow(run))
{
  run[,{
    vl_bsub(cmd,
            cores = 4L,
            m = 20,
            t = "08:00:00",
            name = "vlloub",
            o = "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            e = "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/")
  }, .(cmd, set)]
}