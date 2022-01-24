setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
require(glmnet)
require(randomForest)


#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table.rds")
# Restrict to "smaller" and cleaner twist 12
dat <- dat[library=="T12"]
setkeyv(dat, c("CP", "spacer", "group_L", "group_R"))
# Select dependent & predictor variables
obj <- rbindlist(list(dev= dat[.("DSCP", "SCR1_300")],
                      hk= dat[.("RpS12", "SCR1_300")]), 
                 idcol = "group")
# Add motifs
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
sel_motifs <- unique(feat$motif_padj[padj<1e-5 & log2OR>0, motif])
sel_motifs <- as.character(sel_motifs)
obj <- feat$add_feature(DT= obj, 
                        feature = feat$top_motifs[, c("ID", sel_motifs), with= F])

#-----------------------------------------------#
# MODELLING
#-----------------------------------------------#
# Training variables
vars <- c("median_L", "median_R", 
          paste0(sel_motifs, "_L"), 
          paste0(sel_motifs, "_R"))

#-----------------------------------------------#
# LASSO regression
#-----------------------------------------------#
models <- obj[, {
  # Define train and test subsets
  seed <- 1
  set.seed(seed)
  selL <- sample(unique(L), 0.85*length(unique(L)))
  set.seed(seed)
  selR <- sample(unique(L), 0.85*length(unique(L)))
  train <- L %in% selL & R %in% selR
  # Scaling
  data <- as.matrix(scale(.SD))
  X_train <- data[train,]
  Y_train <- log2FoldChange[train]
  X_test <- data[!train,]
  Y_test <- log2FoldChange[!train]
  
  # Setting alpha = 1 implements lasso regression
  lambdas <- 10^seq(2, -3, by = -.1)
  lasso_reg <- cv.glmnet(X_train,
                         Y_train,
                         alpha = 1,
                         lambda = lambdas,
                         standardize = TRUE,
                         nfolds = 5)
  # Best  lambda
  lambda_best <- lasso_reg$lambda.min
  # Modelling
  model <- glmnet(X_train,
                  Y_train,
                  alpha = 1,
                  lambda = lambda_best,
                  standardize = TRUE)
  # Prediction
  predict_train <- predict(model, 
                           s = lambda_best, 
                           newx = X_train)
  predict_test <- predict(model, 
                          s = lambda_best, 
                          newx = X_test)
  # Evaluation
  eval_train <- vl_model_eval(observed = Y_train,
                              predicted = predict_train)
  setnames(eval_train, function(x) paste0(x, "_train"))
  eval_test <- vl_model_eval(observed = Y_test,
                              predicted = predict_test)
  setnames(eval_test, function(x) paste0(x, "_train"))
  
  # RETURN
  data.table(lasso_model= .(model),
             predict_train= .(predict_train),
             predict_test= .(predict_test),
             eval_train,
             eval_test)
}, group, .SDcols= vars]

