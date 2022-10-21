setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")

#-----------------------------------------------------#
# Linear models on active pairs
#-----------------------------------------------------#
# Define train and test sets
set.seed(1)
dat[dat[, .(set= sample(3)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(3)), R], setR:= i.set, on= "R"]
dat[, set:= .GRP, .(setL, setR)]
# Train linear model for each train set and compute predicted values
model <- lm(formula = log2FoldChange~indL*indR,
            data= dat)
model$CV_rsqs <- dat[, {
  print(set)
  cL <- L
  cR <- R
  train <- dat[!(L %in% cL) & !(R %in% cR)]
  model <- lm(formula = log2FoldChange~indL*indR,
              data= train)
  .(rsq= summary(model)$r.squared)
}, set]

# Compute expected
dat[, additive:= log2(2^indL+2^indR)]
dat[, multiplicative:= indL+indR]
dat[, predicted:= predict(model)]
dat[, residuals:= log2FoldChange-predicted]
dat[, medianResidualsL:= median(residuals), L]
dat[, medianResidualsR:= median(residuals), R]

saveRDS(model, "Rdata/linear_model_vllib002.rds")
