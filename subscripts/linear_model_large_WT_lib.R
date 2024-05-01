setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ---
dat <- readRDS("db/FC_tables/DSCP_large_WT_FC_DESeq2.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Define classes ----
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]

# Train linear model on the whole dataset ----
model <- lm(formula = log2FoldChange ~ indL * indR,
            data= dat)

# Compute expected and residuals ----
dat[, `Linear model`:= predict(model)]
dat[, residuals:= log2FoldChange-`Linear model`]

# Define train and test sets ----
set.seed(1)
dat[dat[, .(set= sample(3)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(3)), R], setR:= i.set, on= "R"]
dat[, set:= .GRP, .(setL, setR)]
dat$setL <- dat$setR <- NULL

# Cross-validated lm ----
model$CV_rsqs <- dat[, {
  print(subset)
  enhL <- L
  enhR <- R
  train <- dat[!(L %in% enhL) & !(R %in% enhR)]
  test <- dat[(L %in% enhL) & (R %in% enhR)]
  model <- lm(formula = log2FoldChange ~ indL * indR,
              data= train)
  # Compute rsq ----
  rsquared <- cor(test$log2FoldChange, predict(model, test))^2
  # Adjusted rsq ----
  nPredictors <- length(coef(model)-1) # Intercept should not be counted
  nObs <- nrow(test)
  .(adj.rsq= 1 - ((1 - rsquared) * (nObs - 1) / (nObs - nPredictors - 1)))
}, set]

# Save ----
saveRDS(model, "db/linear_models/lm_DSCP_large_WT.rds")
saveRDS(dat, "db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")