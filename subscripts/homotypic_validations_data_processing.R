setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(readxl)

# Import and format data ----
dat <- data.table(scheme= list.files("db/luciferase/homotypic_pairs/", "scheme.csv$", full.names = T))
dat[, luc_file:= gsub("_scheme.csv", "_luc.csv", scheme)]
dat[, ren_file:= gsub("_scheme.csv", "_ren.csv", scheme)]
dat <- dat[, {
  .c <- as.data.table(plater::read_plate(scheme))
  .c[, luc:= plater::read_plate(luc_file)$values]
  .c[, ren:= plater::read_plate(ren_file)$values]
  .c$rep <- .GRP
  .c
}, seq(nrow(dat))]
dat[, c("L", "R"):= tstrsplit(rownames, "__")]
dat$seq <- dat$Wells <- dat$rownames <- NULL

# Compute activity ----
dat <- dat[ren>2500]
dat <- dat[, .(log2FoldChange= mean(log2(luc/ren))), .(L, R, rep)]

# Mean Act ----
dat <- dat[, .(log2FoldChange= mean(log2FoldChange),
               sd= sd(log2FoldChange),
               all_vars= .(log2FoldChange)), .(L, R)]

# Center data ----
dat[, log2FoldChange:= log2FoldChange-log2FoldChange[grepl("^control", L) & grepl("^control", R)]]

# Individual act
dat[, indL:= log2FoldChange[grepl("^control", R)], L]
dat[, indR:= log2FoldChange[grepl("^control", L)], R]

# Models ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
dat[, `Linear model`:= predict(model, .SD)]

# Save ----
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
saveRDS(dat, "Rdata/homotypic_validations_luciferase_final_table.rds")