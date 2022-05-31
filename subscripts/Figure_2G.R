setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
screen <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class== "enh./enh."]
screen[, diff:= log2FoldChange-additive]
screen <- feat$add_feature(screen, feat$lib)
cl <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")

# Build data
dat <- screen[L %in% cl$rows$name & R %in% cl$cols$name]

# Sample for CV
set.seed(1)
sel_L <- dat[, 1-.N/dat[,.N], L][, sample(L, round(.N/10), prob = V1)]
set.seed(1)
sel_R <- dat[, 1-.N/dat[,.N], R][, sample(R, round(.N/10), prob = V1)]
dat[, set:= ifelse(L %in% sel_L | R %in% sel_R, "test", "train")]

# Train linear model
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= dat[set=="train"])

# Predict and CV
rsq <- vl_model_eval(observed = dat[set=="test", log2FoldChange], 
                     predicted = predict(model, new= dat[set=="test"]))
pred <- predict(model, new= dat)

#---------------------------------------#
# PLOT
#---------------------------------------#
pdf("pdf/draft/Figure_2G.pdf",
    height = 3,
    width = 3)
par(mgp= c(1.5, 0.5, 0),
    mar= c(3.5,2.5,0.5,1),
    tcl= -0.2,
    yaxs= "r")
smoothScatter(pred, 
              dat$log2FoldChange,
              xlab= "Linear model prediction (log2)",
              ylab= "Observed (log2)",
              colramp = colorRampPalette(c("white", "grey90", "grey60", "grey20")),
              las= 1)
legend("topleft", 
       paste0("CV R²= ", round(rsq$Rsquare, 2)), 
       bty= "n")
abline(0, 1, lty= 2)
dev.off()