setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ---
dat <- readRDS("db/FC_tables/DSCP_ECD_WT_FC_DESeq2.rds")

# Select activity-matched enhancers for OSC and dat ----
ecdL <- vl_select_act_matched(unique(dat[grepl("^ecd", L), .(L, indL)]),
                              unique(dat[grepl("^dev", L), .(L, indL)]),
                              "indL")
ecdL <- ecdL$L
ecdR <- vl_select_act_matched(unique(dat[grepl("^ecd", R), .(R, indR)]),
                              unique(dat[grepl("^dev", R), .(R, indR)]),
                              "indR")
ecdR <- ecdR$R
dat[, group:= fcase(grepl("^ecd", L) & grepl("^ecd", R), "Ecd./Ecd.",
                    grepl("^ecd", L) & R %in% ecdR, "Ecd./S2",
                    L %in% ecdL & grepl("^ecd", R), "S2/Ecd.",
                    L %in% ecdL & R %in% ecdR, "S2/S2")]
dat[, group:= factor(group, c("Ecd./Ecd.", "Ecd./S2", "S2/Ecd.", "S2/S2"))]
dat <- dat[!is.na(group)]

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]

# Train linear model ----
model <- lm(log2FoldChange~0+indL*indR, dat)
dat[, `Linear model`:= predict(model, .SD)]
dat[, residuals:= log2FoldChange-`Linear model`]

saveRDS(model, "db/linear_models/lm_DSCP_ECD_ecd_pair.rds")
saveRDS(dat, "db/linear_models/FC_DSCP_ECD_ecd_pairs_lm_predictions.rds")