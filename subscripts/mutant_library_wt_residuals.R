setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_mutant_library_FC_DESeq2.rds")
dat <- dat[grepl("twist_WT_|Trl_WT_|noMotifAct_WT_", L) & grepl("twist_WT_|Trl_WT_|noMotifAct_WT_", R)]

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("groupL", "mutL", "L"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "R"):= tstrsplit(R, "_", keep= c(1,2,4))]
dat[, L:= paste0(groupL, L)]
dat[, R:= paste0(groupR, R)]

# Define enhancer pairs classes ----
dat[, class:= fcase(grepl("^noMotifAct", L) & grepl("^noMotifAct", R), "No motif",
                    grepl("^twist", L) & grepl("^twist", R), "Twist/Twist",
                    grepl("^Trl", L) & grepl("^Trl", R), "Trl/Trl")]
dat <- dat[!is.na(class)]

# Compute expected multiplicative ----
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
dat[, `linear model`:= predict(model, .SD)]
dat[, residuals:= log2FoldChange-`linear model`]

# Plot
pdf("pdf/draft/residuals_wt_variants.pdf", 6, 3)
vl_par(mai= c(.9, 1.2, .9, 1.2),
       mgp= c(1, .35, 0),
       mfrow= c(1,2))
vl_boxplot(residuals~class,
           dat,
           tilt.names= T,
           compute.pval= list(c(1,2), c(1,3)),
           ylab= "Residuals (log2)")
abline(h= median(dat[class=="No motif", residuals]),
       lty= "11")
vl_boxplot(log2FoldChange~class,
           dat,
           tilt.names= T,
           compute.pval= list(c(1,2), c(1,3)),
           ylab= "Activity (log2)")
abline(h= median(dat[class=="No motif", log2FoldChange]),
       lty= "11")
dev.off()