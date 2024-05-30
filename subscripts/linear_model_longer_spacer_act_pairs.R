setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ---
dat <- readRDS("db/FC_tables/DSCP_long_spacer_FC_DESeq2.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R) & actL!="Inactive" & actR!="Inactive"]

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]

# Train linear model ----
model <- lm(log2FoldChange~0+indL*indR, dat)
dat[, `Linear model`:= predict(model)]
dat[, residuals:= log2FoldChange-`Linear model`]

saveRDS(model, "db/linear_models/lm_DSCP_long_spacer_act_pairs.rds")
saveRDS(dat, "db/linear_models/FC_DSCP_long_spacer_act_pairs_lm_predictions.rds")
