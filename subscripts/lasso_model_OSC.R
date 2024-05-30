setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

dat.file <- "db/linear_models/FC_DSCP_OSC_full_dataset_lm_predictions.rds"

# Import data ---- (test sets already present in the #set column)
dat <- readRDS(dat.file)

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
set[, output:= paste0("db/lasso_models/",
                      ifelse(set==0, "full_dataset", paste0("fold", set)), "_OSC_lasso.rds")]

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
          "git_peSTARRSeq/functions/train_LASSO_model.R",
          excludeL,
          excludeR,
          "db/motif_counts/twist008_motif_counts_selected.rds",
          dat.file,
          output)
  }else
    as.character(NA)
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
            name = "lasOSC",
            o = "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            e = "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/")
  }, .(cmd, set)]
}