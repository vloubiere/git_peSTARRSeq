setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ---
dat <- readRDS("db/FC_tables/DSCP_OSC_WT_FC_DESeq2.rds")

# Select activity-matched enhancers for OSC and dat ----
oscL <- vl_select_act_matched(unique(dat[grepl("^OSC", L), .(L, indL)]),
                              unique(dat[grepl("^dev", L), .(L, indL)]),
                              "indL")
oscL <- oscL$L
oscR <- vl_select_act_matched(unique(dat[grepl("^OSC", R), .(R, indR)]),
                              unique(dat[grepl("^dev", R), .(R, indR)]),
                              "indR")
oscR <- oscR$R
dat[, group:= fcase(grepl("^OSC", L) & grepl("^OSC", R), "OSC/OSC",
                    grepl("^OSC", L) & R %in% oscR, "OSC/S2",
                    L %in% oscL & grepl("^OSC", R), "S2/OSC",
                    L %in% oscL & R %in% oscR, "S2/S2")]
dat[, group:= factor(group, c("OSC/OSC", "OSC/S2", "S2/OSC", "S2/S2"))]
dat <- dat[!is.na(group)]

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]

# Train linear model ----
model <- lm(log2FoldChange~0+indL*indR, dat)
dat[, `Linear model`:= predict(model)]
dat[, residuals:= log2FoldChange-`Linear model`]

saveRDS(model, "db/linear_models/lm_DSCP_OSC_WT.rds")
saveRDS(dat, "db/linear_models/FC_DSCP_OSC_WT_lm_predictions.rds")