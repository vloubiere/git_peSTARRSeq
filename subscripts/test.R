setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
old <- readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds")
dat[old, log2FoldChange.old:= i.log2FoldChange, on= c("L", "R")]
setnames(dat, c("mean_L", "mean_R"), c("indL", "indR"))
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
dat[, nameL:= LETTERS[.GRP], L]
dat[, nameR:= LETTERS[.GRP], R]
dat[, additive:= log2(2^indL+2^indR-1)]
# model <- lm(log2FoldChange~indL*indR, dat)
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
summary(model)
# dat[, multiplicative:= predict(model)]
dat[, multiplicative:= predict(model, .SD)]
# dat[, multiplicative:= indL+indR]
setorderv(dat, "log2FoldChange")
barplot(t(dat[, .(additive, log2FoldChange, multiplicative)]),
        beside = T,
        col= c("grey80", "grey10", "tomato"),
        border= NA)

barplot(dat$additive, space= 2, col= "lightgrey", border= NA, )
barplot(dat$log2FoldChange, space= 2, col= "tomato", border= NA, add= T)

model <- lm(log2FoldChange~I(2^indL)*I(2^indR), dat)
# model <- lm(log2FoldChange~indL*indR, dat)
summary(model)
dat[, multiplicative:= predict(model)]
dat[, best:= ifelse(abs(log2FoldChange-additive)<abs(log2FoldChange-multiplicative), "Additive", "Multiplicative")]

par(mfrow= c(1,2))
plot(dat$additive, dat$log2FoldChange)
abline(0,1)
plot(dat$multiplicative, dat$log2FoldChange)
abline(0,1)

