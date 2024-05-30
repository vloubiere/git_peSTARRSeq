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
require(glmnet)
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

# Counts matrix (sum left and right counts) ----
motL <- counts[data$L, on= "ID"][, -1]
motR <- counts[data$R, on= "ID"][, -1]
counts <- as.matrix(motL)+as.matrix(motR)

# CLean memory ----
rm(list= c("motL", "motR"))
gc()

# Train model ----
res <- list()
for(response in c("log2FoldChange", "residuals"))
{
  # Setting alpha = 1 implements lasso regression ----
  lambdas <- 10^seq(2, -3, by = -.1)
  lasso_reg <- cv.glmnet(counts[data$set=="train",],
                         data[set=="train"][[response]],
                         alpha = 1,
                         lambda = lambdas,
                         standardize = TRUE,
                         nfolds = 5)
  # Best  lambda ----
  lambda_best <- lasso_reg$lambda.min
  
  # Modelling ----
  model <- glmnet(counts[data$set=="train",],
                  data[set=="train"][[response]],
                  alpha = 1,
                  lambda = lambda_best,
                  standardize = TRUE)
  
  # Predict ----
  predVar <- paste0("lasso_predicted_", response)
  data[, (predVar):= predict(model,
                             s= model$lambda_best,
                             newx= counts)]
  
  # Compute rsq ----
  rsquared <- if(any(data$set=="test"))
    cor(data[set=="test"][[response]], data[set=="test"][[predVar]])^2 else 
      cor(data[[response]], data[[predVar]])^2# If not test set, rsq will be quantified on the whole dataset
  
  # Adjusted rsq ----
  coefficients <- coef(model, s = model$best_lambda)
  nPredictors <- sum(coefficients != 0) - 1 # Intercept should not be counted
  nObs <- if(any(data$set=="test")) sum(data$set=="test") else nrow(data)
  model$adj.rsq <- 1 - ((1 - rsquared) * (nObs - 1) / (nObs - nPredictors - 1))
  
  # Save ----
  res[[response]] <- model
}
res[["data"]] <- data

# Save results ----
saveRDS(res, output)