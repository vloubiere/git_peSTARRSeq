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
dat <- dat[vllib %in% c("vllib016", "vllib015")][group_L!="control" & group_R!="control"]
# dat <- dat[vllib=="vllib016"][group_L!="control" & group_R!="control"]

# subgroups
set.seed(1)
dat[, subset_L:= sample(c("train", "test"), .N, replace = T, prob = c(1.5, 1)), L]
set.seed(2)
dat[, subset_R:= sample(c("train", "test"), .N, replace = T, prob = c(1.5, 1)), R]
length(which(dat$subset_L=="test" & dat$subset_R=="test"))/nrow(dat)
dat[, subset:= ifelse(subset_L=="test" & subset_R=="test", "test", "train")]
dat[, col:= colorRampPalette(unlist(.BY))(3)[2], .(col_L, col_R)]
setkeyv(dat, "subset")

data <- dat[, .(CP, subset, log2FoldChange, sapply(.SD, scale)), .SDcols= c("median_L", "median_R")]
.lm <- lm(log2FoldChange~CP+median_L*median_R,
          data["train"])

# Prediction
predict_train <- predict(.lm)
predict_test <- predict(.lm,
                        newdata = data["test"])
# Evaluation
eval_train <- vl_model_eval(observed = data["train", log2FoldChange],
                            predicted = predict_train)
setnames(eval_train, function(x) paste0(x, "_train"))
eval_test <- vl_model_eval(observed = data["test", log2FoldChange],
                           predicted = predict_test)
setnames(eval_test, function(x) paste0(x, "_test"))
eval_train
eval_test


obs <- data$log2FoldChange
exp <- predict(.lm, newdata= data)

# pdf("test/test_scatterplots.pdf")
par(mfrow= c(2,3),
    las= 1,
    pch= 16)
plot(obs,
     exp,
     col= adjustcolor(dat$col, 0.5),
     cex= 0.5)
abline(0, 1)
legend("topleft", 
       paste0("PCC= ", 
              round(cor.test(obs, exp)$estimate, 2), 
              "\nR2= ", 
              round(vl_model_eval(obs,
                                  exp)$Rsquare, 2)),
       bty= "n")
legend("bottomright", 
       paste0("PCC= ", 
              round(cor.test(obs, exp)$estimate, 2), 
              "\nR2= ", 
              round(vl_model_eval(obs,
                                  exp)$Rsquare, 2)),
       bty= "n")




plot(data$log2FoldChange,
     data$log2FoldChange-predict(.lm, newdata= data),
     col= adjustcolor(dat$col, 0.5),
     pch= 16,
     cex= 0.7)
abline(h= 0)
dev.off()




setkeyv(dat, c("CP", "spacer", "group_L", "group_R"))
# Select dependent & predictor variables
obj <- rbindlist(list(dev= dat[.("DSCP", "SCR1_300")],
                      hk= dat[.("RpS12", "SCR1_300")]), 
                 idcol = "group")
# Add motifs
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
sel_motifs <- unique(feat$motif_padj[padj<1e-5 & log2OR>0, motif])
sel_motifs <- as.character(sel_motifs)
new_sel_motifs <- rename_motif_cols(sel_motifs)
setnames(feat$top_motifs, 
         sel_motifs, 
         new_sel_motifs)
obj <- feat$add_feature(DT= obj, 
                        feature = feat$top_motifs[, c("ID", new_sel_motifs), with= F])









#---------------------------------------------#
# Modelling
#---------------------------------------------#
predictors <- readRDS("Rdata/lasso_predictors_selection.rds")
predictors <- predictors$top_15_predictors[!best_predictors %in% c("median_L", "median_R", "(Intercept)")]
predictors[, best_predictors:= rename_motif_cols(best_predictors)]
setkeyv(predictors, "group")

# subgroups
set.seed(1)
obj[, subset_L:= sample(c("train", "test"), .N, replace = T, prob = c(1.5, 1)), L]
set.seed(2)
obj[, subset_R:= sample(c("train", "test"), .N, replace = T, prob = c(1.5, 1)), R]
length(which(obj$subset_L=="test" & obj$subset_R=="test"))/nrow(obj)
obj[, subset:= ifelse(subset_L=="test" & subset_R=="test", "test", "train")]

# modelling
models <- obj[, {
  # Subsets
  vars <- predictors[group, best_predictors]
  data <- .SD[, .(log2FoldChange, sapply(.SD, scale)), .SDcols= c("median_L", "median_R", vars)]
  train <- data[subset=="train"]
  test <- data[subset=="test"]
  
  # Modelling
  model <- lm(paste0("log2FoldChange~median_L*median_R+",
                     paste0(vars, collapse= "+")),
              data = train)
  # Prediction
  predict_train <- predict(model)
  predict_test <- predict(model,
                          newdata = test)
  predict_all <- predict(model,
                         newdata = data)
  # Evaluation
  eval_train <- vl_model_eval(observed = train$log2FoldChange,
                              predicted = predict_train)
  setnames(eval_train, function(x) paste0(x, "_train"))
  eval_test <- vl_model_eval(observed = test$log2FoldChange,
                             predicted = predict_test)
  setnames(eval_test, function(x) paste0(x, "_test"))

  # RETURN
  data.table(model= .(model),
             train= .(train),
             test= .(test),
             predict_train= .(predict_train),
             predict_test= .(predict_test),
             eval_train,
             eval_test,
             log2FoldChange= .(data$log2FoldChange),
             predict_all= .(predict_all))
}, group]








