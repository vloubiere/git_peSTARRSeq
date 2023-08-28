setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/vllib029_DESeq2.rds")
dat[, c("groupL", "mutL", "IDL"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "IDR"):= tstrsplit(R, "_", keep= c(1,2,4))]
# Make sequence ID unique (numbers are only unique within one group)
dat[, IDL:= paste0(groupL, IDL)]
dat[, IDR:= paste0(groupR, IDR)]
dat[, multiplicative:= indL+indR]
dat[, residuals:= log2FoldChange-multiplicative]

saveRDS(dat,
        "db/linear_models/FC_vllib029_mutLib_lm_predictions.rds")