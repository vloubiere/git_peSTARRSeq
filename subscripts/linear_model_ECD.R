setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ---
dat <- readRDS("db/FC_tables/DSCP_ECD_WT_DESeq2.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Linear models on all pairs ----
## Define train and test sets ----
set.seed(1)
dat[dat[, .(set= sample(3)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(3)), R], setR:= i.set, on= "R"]
dat[, set:= .GRP, .(setL, setR)]
## Train linear model for each train set and compute predicted values ----
model <- lm(formula = log2FoldChange ~ indL * indR,
            data= dat)
model$CV_rsqs <- dat[, {
  print(set)
  cL <- L
  cR <- R
  train <- dat[!(L %in% cL) & !(R %in% cR)]
  model <- lm(formula = log2FoldChange ~ indL * indR,
              data= train)
  .(rsq= summary(model)$r.squared)
}, set]

## Compute expected ----
dat[, predicted:= predict(model)]
dat[, residuals:= log2FoldChange-predicted]

saveRDS(model, "db/linear_models/lm_DSCP_ECD_large_WT.rds")
saveRDS(dat, "db/linear_models/FC_DSCP_ECD_large_WT_lm_predictions.rds")