setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data and remove control pairs ----
dat <- readRDS("db/FC_tables/DSCP_ECD_WT_FC_DESeq2.rds")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Select activity-matched enhancers for ECD ----
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

# Predicted values ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]
model <- lm(log2FoldChange ~ indL * indR, dat)
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
        "db/linear_models/lm_DSCP_ECD_full_dataset.rds")
saveRDS(dat,
        "db/linear_models/FC_DSCP_ECD_full_dataset_lm_predictions.rds")