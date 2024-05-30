setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)

# Import data and remove control pairs ----
dat <- readRDS("db/FC_tables/DSCP_OSC_WT_FC_DESeq2.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Predicted values ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]
model <- lm(log2FoldChange~indL*indR, dat)
dat[, `Linear model`:= predict(model)]
dat[, residuals:= log2FoldChange-`Linear model`]

# Define train and test sets ----
set.seed(1)
dat[dat[, .(set= sample(3)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(3)), R], setR:= i.set, on= "R"]
dat[, set:= .GRP, .(setL, setR)]
dat$setL <- dat$setR <- NULL

saveRDS(model,
        "db/linear_models/lm_DSCP_OSC_full_dataset.rds")
saveRDS(dat,
        "db/linear_models/FC_DSCP_OSC_full_dataset_lm_predictions.rds")