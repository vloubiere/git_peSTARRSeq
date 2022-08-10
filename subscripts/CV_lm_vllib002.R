setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act== "enh./enh."]
set.seed(1)
dat <- lib[sample(nrow(lib), nrow(lib))]

# Define train and test sets
set.seed(1)
sel_L <- dat[, 1-.N/dat[,.N], L][, sample(L, round(.N/5), prob = V1)]
set.seed(1)
sel_R <- dat[, 1-.N/dat[,.N], R][, sample(R, round(.N/5), prob = V1)]
dat[L %in% sel_L & R %in% sel_R, set:= "test"]
dat[!(L %in% sel_L) & !(R %in% sel_R), set:= "train"]

# Train linear model
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= dat[set=="train"])
dat[, predicted:= predict(model, 
                          newdata = dat)]
rsq <- vl_model_eval(observed = dat[set=="test", log2FoldChange], 
                     predicted = predict(model, new= dat[set=="test"]))
eq <- vl_model_equation(model, digits = 2)
eq <- gsub("median_L", "5'", eq)
eq <- gsub("median_R", "3'", eq)
eq <- gsub("log2FoldChange ", "Act.", eq)
eq <- gsub("\\s\\*\\s", "*", eq)
saveRDS(list(vllib002_dat= dat,
             model= model,
             rsq= rsq,
             eq= eq), 
        "Rdata/CV_linear_model_vllib002.rds")