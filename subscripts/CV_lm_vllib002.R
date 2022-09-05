setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
# lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act== "enh./enh."]
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
set.seed(1)
dat <- lib[sample(nrow(lib), nrow(lib))]

# Define train and test sets
set.seed(1)
dat[dat[, .(set= sample(5)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(5)), R], setR:= i.set, on= "R"]
dat[, set:= paste0("test", .GRP), .(setL, setR)]
# Train linear model for each train set and compute predicted values
dat[, c("predicted", "intercept", "coeff_L", "coeff_R", "coeff_LR", "R2"):= {
  print(set)
  cL <- L
  cR <- R
  train <- dat[!(L %in% cL) & !(R %in% cR)]
  model <- lm(formula = log2FoldChange~median_L*median_R,
              data= train)
  pred <- predict(model, newdata = .SD)
  rsq <- vl_model_eval(observed = log2FoldChange, 
                       predicted = pred)$Rsquare
  print(vl_model_equation(model))
  c(list(pred), as.list(model$coefficients), list(rsq))
}, set]

model <- lapply(unique(dat[, .(intercept, coeff_L, coeff_R, coeff_LR, R2)]), mean)
model <- list(equation= paste0("Activity= ", 
                               round(model$intercept, 2),
                               "+",
                               round(model$coeff_L, 2), "*5'+",
                               round(model$coeff_R, 2), "*3'",
                               round(model$coeff_LR, 2), "*5':3'"),
              R2= model$R2)
saveRDS(list(pred= dat[, .(L, R, predicted, R2)],
             model= model), 
        "Rdata/CV_linear_model_vllib002.rds")
