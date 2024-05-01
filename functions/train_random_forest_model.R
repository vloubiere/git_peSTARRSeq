#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############--------------------------------------------------------------------##############
############  Compute individual and paired activities from pSTARR-Seq counts   ##############
############--------------------------------------------------------------------##############
# Takes as input umi collapsed counts and computes log2FoldChange, using either DESeq2 approach or log2 ratios

# test if there is at least 2 args: if not, return an error
if (length(args)!=5) {
  stop(paste0(length(args), " / 5 required arguments should be provided:\n",
              paste0(args, "\n"),
              "Please specify:\n
              [required] 1/ Comma-separated list of left enhancers to exclude from training\n
              [required] 2/ Comma-separated list of right enhancers to exclude from training\n
              [required] 3/ Motif counts matrix with first column containing enhancer IDs \n
              [required] 4/ .rds input data file containing the response variables \n
              [required] 5/ .rds output file \n"))
}
require(randomF)
require(data.table)
require(vlfunctions)

# Test args ----
data <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
excludeL <- sample(unique(data$L), size = 100)
excludeR <- sample(unique(data$R), size = 100)
counts <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/db/motif_counts/twist008_motif_counts_selected.rds")

# # Extract variables ----
# excludeL <- unique(unlist(tstrsplit(args[1], ",")))
# excludeR <- unique(unlist(tstrsplit(args[2], ",")))
# counts <- readRDS(args[3])
# data <- readRDS(args[4])
# output <- args[5]

# Split training and test datasets ----
data[, set:= fcase(!(L %in% excludeL) & !(R %in% excludeR), "train",
                   (L %in% excludeL) & (R %in% excludeR), "test",
                   default = "notUsed")]

# Counts matrix ----
motL <- counts[data$L, on= "ID"][, -1]
setnames(motL, function(x) paste0(x, "__L"))
motR <- counts[data$R, on= "ID"][, -1]
setnames(motR, function(x) paste0(x, "__R"))
counts <- cbind(motL, motR)
counts <- as.matrix(counts)
rm(list= c("motL", "motR"))
gc()

# Train model ----
res <- list()
for(response in c("log2FoldChange", "residuals"))
{
  model <- randomForest::randomForest(x= counts[data$set=="train",],
                                      y= data[set=="train"][[response]],
                                      importance= TRUE,
                                      proximity= TRUE)
  
  # Predict ----
  predVar <- paste0("rf_predicted_", response)
  data[, (predVar):= predict(model,
                             s= model$lambda_best,
                             newx= counts)]
  
  # Save ----
  model$rsq <- if(any(data$set=="test"))
    vl_model_eval(data[set=="test"][[response]],
                  data[set=="test"][[predVar]])$Rsquare else
                    vl_model_eval(data[[response]], data[[predVar]])$Rsquare
  res[[response]] <- model
}
res[["data"]] <- data

# Save results ----
saveRDS(res,
        output)