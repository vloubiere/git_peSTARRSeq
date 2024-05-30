setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ---
dat <- readRDS("db/FC_tables/DSCP_focused_WT_FC_DESeq2.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Define classes ----
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]

# Save ----
saveRDS(dat,
        "db/linear_models/FC_DSCP_focused_lm_predictions.rds")