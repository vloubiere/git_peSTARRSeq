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
require(randomForest)
require(data.table)
require(vlfunctions)

# Extract variables ----
excludeL <- unique(unlist(tstrsplit(args[1], ",")))
excludeR <- unique(unlist(tstrsplit(args[2], ",")))
counts <- readRDS(args[3])
data <- readRDS(args[4])
output <- args[5]

# Split training and test datasets ----
data[, set:= fcase(!(L %in% excludeL) & !(R %in% excludeR), "train",
                   (L %in% excludeL) & (R %in% excludeR), "test",
                   default = "notUsed")]
# data <- data[set=="train"][1:200] # For tests

# Counts matrix (keep left and right counts separated) ----
motL <- counts[data$L, on= "ID"][, -1]
setnames(motL, function(x) paste0(x, "__L"))
motR <- counts[data$R, on= "ID"][, -1]
setnames(motR, function(x) paste0(x, "__R"))
counts <- cbind(motL, motR)
# Sum values left and right counts)
# motL <- counts[data$L, on= "ID"][, -1]
# motR <- counts[data$R, on= "ID"][, -1]
# counts <- as.matrix(motL)+as.matrix(motR)

# CLean memory ----
rm(list= c("motL", "motR"))
gc()

# Train model ----
res <- list()
for(response in c("log2FoldChange", "residuals"))
{
  # Train model ----
  model <- randomForest(x= counts[data$set=="train",],
                        y= data[set=="train"][[response]],
                        importance= TRUE)
  
  # Predict ----
  predVar <- paste0("rf_predicted_", response)
  data[, (predVar):= predict(model, newdata= counts)]
  
  # Save ----
  res[[response]] <- model
}
res[["data"]] <- data

# Save results ----
saveRDS(res, output)