setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_mutant_library_FC_DESeq2.rds")

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("groupL", "mutL", "L"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "R"):= tstrsplit(R, "_", keep= c(1,2,4))]
dat[, L:= paste0(groupL, L)]
dat[, R:= paste0(groupR, R)]
dat <- dat[!(grepl("^control", L) & mutL=="WT") & !(grepl("^control", R) & mutR=="WT")]

# Compute expected multiplicative ----
model <- lm(log2FoldChange~indL*indR, dat[grepl("WT", mutL) & grepl("WT", mutR)])
dat[, `linear model`:= predict(model, .SD)]
dat[, `additive model`:= log2(2^indL+2^indR-1)]
# dat[, residuals:= log2FoldChange-`additive model`]
dat[, residuals:= log2FoldChange-`linear model`]
dat <- dat[padjL<0.05 & indL>log2(1.5) & padjR<0.05 & indR>log2(1.5)]

dat[,cdition:= fcase(mutL=="WT" & mutR=="WT", "WT",
                     grepl("add.*Trl", mutL) & grepl("add.*Trl", mutR) & grepl("^noMotifAct", L) & grepl("^noMotifAct", R), "Added Trl motifs",
                     grepl("add.*Twist", mutL) & grepl("add.*Twist", mutR) & grepl("^noMotifAct", L) & grepl("^noMotifAct", R), "Added Twist motifs",
                     grepl("add.*Dref", mutL) & grepl("add.*Dref", mutR), "Added Dref motifs",
                     grepl("mut.*Trl", mutL) & grepl("mut.*Trl", mutR), "Mutated Trl motifs",
                     grepl("mut.*Twist", mutL) & grepl("mut.*Twist", mutR), "Mutated Twist motifs")]
dat <- dat[!is.na(cdition)]

test <- merge(dat[cdition=="WT", !"cdition"],
              dat[cdition!="WT"],
              by= c("L", "R"),
              suffixes= c("_wt", "_mut"))
test <- test[abs(log2FoldChange_wt-log2FoldChange_mut)<.5]
test <- melt(test,
             id.vars = "cdition",
             measure.vars = c("residuals_wt", "residuals_mut"))
vl_boxplot(value~variable+cdition,
           test,
           names= function(x) gsub("residuals", "", x),
           tilt.names= T,
           compute.pval= list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10)),
           notch= T)
