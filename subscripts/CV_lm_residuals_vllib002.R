setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh."]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
lib[, predicted:= predict(model, newdata = .SD)]
lib[, residuals:= log2FoldChange-predicted]

# Train models to predict left and right median residuals
if(!file.exists("Rdata/CV_LASSO_residuals_vllib002.rds"))
{
  models <- list()
  feat <- fread("Rdata/final_300bp_enhancer_features.txt")
  for(pos in c("L", "R"))
  {
    # Retrieve data
    .c <- lib[, .(residuals= median(residuals)), c(pos, paste0("median_", pos))]
    setnames(.c, c("ID", "median", "residuals"))
    .c <- cbind(.c, feat[.c, on= "ID"])
    
    # Split predictors and dependent var
    x <- cbind(.c[, .(median)],
               .c[, names(.c) %in% vl_Dmel_motifs_DB_full$motif, with= F])
    y <- .c$residuals
    
    ########################
    # Full model
    ########################
    # perform k-fold cross-validation to find optimal lambda value
    set.seed(1)
    cv_model <- cv.glmnet(as.matrix(x), y, alpha = 1)
    # find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    # Model
    set.seed(1)
    best_model <- glmnet(x= as.matrix(x), 
                         y= y, 
                         alpha = 1, 
                         lambda = best_lambda)
    # Assess performance
    predicted <- predict(best_model, newx = as.matrix(x))
    # Return
    rsq <- vl_model_eval(y, predicted)$Rsquare
    
    ########################
    # CV sub models
    ########################
    CV_rsq <- mclapply(1:25, function(i)
    {
      set.seed(i)
      sample <- sample(length(y), round(length(y)*0.9))
      cv_model <- cv.glmnet(as.matrix(x)[sample,], y[sample], alpha = 1)
      best_lambda <- cv_model$lambda.min
      best_model <- glmnet(x= as.matrix(x)[sample,], 
                           y= y[sample], 
                           alpha = 1, 
                           lambda = best_lambda)
      predicted <- predict(best_model, newx = as.matrix(x)[-sample,])
      print(i)
      # Return
      vl_model_eval(y[-sample], predicted)$Rsquare
    }, mc.preschedule = F, mc.cores = getDTthreads())
    
    
    models[[pos]] <- list(model= best_model,
                          coeffs= as.data.table(as.matrix(coef(best_model)), 
                                                keep.rownames = var),
                          rsq= rsq,
                          CV_rsq= unlist(CV_rsq))
  }
  
  saveRDS(models,
          "Rdata/CV_LASSO_residuals_vllib002.rds")
}

smoothScatter(models$,
              lib$residuals)

lib[, strength_L:= cut(median_L, c(-Inf, 2,4,6, Inf), c("Weak", "Med.", "Strong", "Very Str."))]
lib[, strength_R:= cut(median_R, c(-Inf, 2,4,6, Inf), c("Weak", "Med.", "Strong", "Very Str."))]
mat <- dcast(lib, strength_L~strength_R, value.var = "residuals", fun.aggregate = mean)
mat <- mat[nrow(mat):1,]
vl_heatmap(as.matrix(mat, 1), 
           cluster_rows = F, 
           cluster_cols = F,
           breaks= c(-1,0,1))

lib[feat, DRE_L:= `DRE/1__2`, on= "L==ID"]
lib[DRE_L>5, DRE_L:= 5]
lib[feat, DRE_R:= `DRE/1__2`, on= "R==ID"]
lib[DRE_R>5, DRE_R:= 5]
mat <- dcast(lib, DRE_L~DRE_R, value.var = "residuals", fun.aggregate = mean)
mat <- mat[nrow(mat):1,]
vl_heatmap(as.matrix(mat, 1), 
           cluster_rows = F, 
           cluster_cols = F,
           breaks= c(-1,0,1))

lib[feat, AP1_L:= `HD/16__6`, on= "L==ID"]
lib[AP1_L>3, AP1_L:= 3]
lib[feat, AP1_R:= `HD/16__6`, on= "R==ID"]
lib[AP1_R>3, AP1_R:= 3]
mat <- dcast(lib, AP1_L~AP1_R, value.var = "residuals", fun.aggregate = mean)
mat <- mat[nrow(mat):1,]
mat <- as.matrix(mat, 1)
mat[is.nan(mat)] <- NA
vl_heatmap(mat, 
           cluster_rows = F, 
           cluster_cols = F,
           breaks= c(-1,0,1))




