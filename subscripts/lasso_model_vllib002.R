setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)

# Import data ---- (test sets already present in the #set column)
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")

# Counts matrix ----
count <- readRDS("db/motif_counts/twist008_motif_counts.rds")
motL <- count[dat$L, on= "ID"]
setnames(motL, gsub("$", "__L", names(motL)))
motR <- count[dat$R, on= "ID"]
setnames(motR, gsub("$", "__R", names(motR)))
mat <- cbind(as.matrix(motL, 1),
             as.matrix(motR, 1))
rm(list= c("motL", "motR"))
gc()

# CROSS VALIDATED LASSO  ----
# Function
trainLASSO <- function(data, 
                       response= var,
                       counts)
{
  # Setting alpha = 1 implements lasso regression
  lambdas <- 10^seq(2, -3, by = -.1)
  lasso_reg <- cv.glmnet(counts,
                         data[[response]],
                         alpha = 1,
                         lambda = lambdas,
                         standardize = TRUE,
                         nfolds = 5)
  # Best  lambda
  lambda_best <- lasso_reg$lambda.min
  # Modelling
  model <- glmnet(counts,
                  data[[response]],
                  alpha = 1,
                  lambda = lambda_best,
                  standardize = TRUE)
  return(model)
}

# Train model ---- with the full dataset and compute expected
models <- list()
for(var in c("log2FoldChange", "residuals"))
{
  model <- trainLASSO(data= dat,
                      response = var,
                      counts= mat)
  predVar <- switch(var, 
                    "log2FoldChange"= "predicted_lasso",
                    "residuals"= "predicted_residuals_lasso")
  dat[, (predVar):= predict(model,
                            s= model$lambda_best,
                            newx= mat)]
  # Train model for each train set and evaluate model (CV)
  model$CV_rsqs <- dat[, {
    print(set)
    # Subset data
    train <- !(dat$L %in% L) & !(dat$R %in% R)
    test <- dat$L %in% L & dat$R %in% R
    # Model
    model <- trainLASSO(data= dat[(train)],
                        response = var,
                        counts= mat[train, ])
    # Evaluation
    predict_test <- predict(model, 
                            s = model$lambda_best, 
                            newx = mat[test, ])
    eval_test <- vl_model_eval(observed = dat[(test)][[var]],
                               predicted = predict_test)
    .(rsq= eval_test$Rsquare)
  }, set]
  models[[var]] <- model
}

saveRDS(models, "db/linear_models/lasso_vllib002_residuals.rds")
saveRDS(dat, "db/linear_models/FC_lasso_vllib002_residuals_predictions.rds")
