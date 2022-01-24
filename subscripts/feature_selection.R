setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
require(glmnet)

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
# Add chromatin features and replace NAs with 0
obj <- feat$add_feature(DT= obj, 
                        feature = feat$chromatin_features)

# Training variables
vars <- c("median_L", "median_R", 
          paste0(sel_motifs, "_L"), 
          paste0(sel_motifs, "_R"))

#-----------------------------------------------#
# LASSO regression
#-----------------------------------------------#
# scaling
obj[, (vars):= lapply(.SD, scale), .SDcols= vars]
# subgroups
set.seed(1)
obj[, subset_L:= sample(seq(4), .N, replace = T), L]
set.seed(2)
obj[, subset_R:= sample(seq(4), .N, replace = T), R]
obj[, subset:= .GRP, .(subset_L, subset_R)]
# modelling
models <- obj[subset!=16, {
  res <- list()
  for(i in 1:15)
  {
    # Subsets
    X_train <- as.matrix(.SD[subset!=i,])
    Y_train <- log2FoldChange[subset!=i]
    X_test <- as.matrix(.SD[subset==i,])
    Y_test <- log2FoldChange[subset==i]
    
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
    setnames(eval_test, function(x) paste0(x, "_test"))
    
    # RETURN
    res[[i]] <- data.table(lasso_model= .(model),
                           predict_train= .(predict_train),
                           predict_test= .(predict_test),
                           eval_train,
                           eval_test)
    print("DONE")
  }
  rbindlist(res)
}, group, .SDcols= vars]

pdf("test/best_predictors.pdf", width = 7)
par(mfrow= c(2,2),
    mar= c(8,4,2,2))
select <- models[, {
  predictors <- rbindlist(lapply(lasso_model, function(x) as.data.table(as.matrix(coef(x)), keep.rownames = T)))
  predictors <- predictors[, .(mean_coeff= mean(s0)), rn]
  predictors[, order:= abs(mean_coeff)]
  setorderv(predictors, "order", -1)
  predictors[15:.N, Cc:= "grey"]
  predictors[1:15, Cc:= ifelse(mean_coeff>0, "tomato", "cornflowerblue")]
  setorderv(predictors, "mean_coeff")
  
  barplot(predictors$mean_coeff,
          las= 2,
          border= NA,
          col= predictors$Cc)
  predictors <- predictors[Cc!="grey"]
  barplot(predictors$mean_coeff,
          las= 2,
          border= NA,
          col= predictors$Cc, 
          names.arg= predictors$rn,
          las= 2)
  .(best_predictors= predictors$rn)
}, group]
dev.off()

obj <- list(lasso_modelling= models,
            top_15_predictors= select)

saveRDS(obj, "Rdata/lasso_predictors_selection.rds")