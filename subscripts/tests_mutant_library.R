setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_mutant_library_FC_DESeq2.rds")
dat <- dat[, -c("groupL", "groupR", "padj")]

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("groupL", "mutL", "IDL"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "IDR"):= tstrsplit(R, "_", keep= c(1,2,4))]
dat[, L:= paste0(groupL, IDL)]
dat[, R:= paste0(groupR, IDR)]
dat$IDL <- dat$IDR <- NULL

# Compute additive and multiplicative expected scores ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]
# Linear model ----
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
dat[, `Linear model`:= predict(model, .SD)]
dat[, residuals:= log2FoldChange-predict(model, .SD)]

# Extract data for combinations of interest ----
dat <- merge(dat[mutL!="WT" | mutR!="WT"],
            dat[mutL=="WT" & mutR=="WT", !c("mutL", "mutR")],
            by= c("L", "R"),
            suffixes= c("", "_wt"))
dat[, mutL:= fcase(mutL=="add2Dref", "addDref",
                  mutL %in% c("add2Trl", "add3Trl"), "addTrl",
                  mutL %in% c("add2Twist", "add3Twist"), "addTwist",
                  default = mutL), mutL]
dat[, mutR:= fcase(mutR=="add2Dref", "addDref",
                  mutR %in% c("add2Trl", "add3Trl"), "addTrl",
                  mutR %in% c("add2Twist", "add3Twist"), "addTwist",
                  default = mutR), mutR]

# Select meaningful combinations ----
dat[grepl("^noMotifAct", L) & grepl("^noMotifAct", R),
    mutClass:= fcase(mutL=="addTrl"   & mutR=="addTwist", "Add Trl/Twist",
                  mutL=="addTwist" & mutR=="addTrl",   "Add Twist/Trl",
                  
                  mutL=="addTwist" & mutR=="addTwist", "Add Twist/Twist",
                  mutL=="addTwist" & mutR=="WT",       "Add Twist/wt",
                  mutL=="WT"       & mutR=="addTwist", "Wt/add Twist",
                  
                  mutL=="addTrl"   & mutR=="addTrl",   "Add Trl/Trl",
                  mutL=="addTrl"   & mutR=="WT",       "Add Trl/wt",
                  mutL=="WT"       & mutR=="addTrl",   "Wt/add Trl",
                  
                  mutL=="addDref"   & mutR=="addDref", "Add Dref/Dref",
                  mutL=="addDref"   & mutR=="WT",      "Add Dref/wt",
                  mutL=="WT"       & mutR=="addDref",  "Wt/add Dref")]
dat[grepl("^mut", L) & grepl("^mut", R),
    mutClass:= fcase(mutL=="mutTrl"   & mutR=="mutTwist", "Trl/Twist mutant",
                  mutL=="mutTwist" & mutR=="mutTrl",   "Twist/Trl mutant",
                  
                  mutL=="mutTwist" & mutR=="mutTwist", "Twist/Twist mutant",
                  mutL=="mutTrl"   & mutR=="mutTrl",   "Trl/Trl mutant",
                  
                  mutL=="mutDref"  & mutR=="mutDref",  "Dref/Dref mutant")]
dat[, wtClass:= fcase(grepl("noMotifAct", L) & grepl("noMotifAct", R), "No motif/No motif",
                      grepl("Trl", L) & grepl("Trl", R), "Trl/Trl",
                      grepl("twist", L) & grepl("twist", R), "Twist/Twist",
                      grepl("noMotifAct", L) & grepl("Trl", R), "No motif/Trl",
                      grepl("noMotifAct", L) & grepl("twist", R), "No motif/Twist",
                      grepl("Trl", L) & grepl("noMotifAct", R), "Trl/No motif",
                      grepl("twist", L) & grepl("noMotifAct", R), "Twist/No motif",
                      grepl("Trl", L) & grepl("twist", R), "Trl/Twist",
                      grepl("twist", L) & grepl("Trl", R), "Twist/Trl")]
dat[, actWtPair:= padjL_wt<0.05 & indL_wt>1 & padjR_wt<0.05 & indR_wt>1]
dat[, actMutPair:= padjL<0.05 & indL>1 & padjR<0.05 & indR>1]

pdf("pdf/draft/review_delta_FC_vs_delta_residuals.pdf", 6, 3)
vl_par(mfrow= c(1,2))
dat[!is.na(mutClass) & (actMutPair), {
  plot(`Linear model`-`Linear model_wt`,
       log2FoldChange-log2FoldChange_wt,
       pch= 16,
       # col= adjustcolor(col, .5),
       cex= 0.25,
       main= mutClass,
       xlab= "Mutant / wt residuals (log2)",
       ylab= "Mutant / wt activity (log2)")
  abline(0, 1, lty= "11")
  plot(`Linear model`,
       log2FoldChange,
       pch= 16,
       # col= adjustcolor(col, .5),
       cex= 0.25,
       main= mutClass)
  abline(0, 1, lty= "11")
}, mutClass]
dev.off()
