setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act== "enh./enh."]
set.seed(1)
dat <- lib[sample(nrow(lib), nrow(lib))]
# Add residuals
res <- readRDS("Rdata/CV_linear_model_vllib002.rds")$pred
dat[res, residuals:= log2FoldChange-i.predicted, on= c("L", "R")]
# Select motifs that were enriched between clusters
enr <- readRDS("Rdata/vllib002_lm_residuals_SOM_motifs_enr.rds")
enr <- rbindlist(enr)[padj<0.00001]
setorderv(enr, "padj")
sel <- unique(enr[rowid(enr$cl)<=10, motif_ID])
# Add motifs to table
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
feat <- feat[, c("ID", paste0(sel, "_motif")), with= F]
setnames(feat, gsub("_motif$", "", names(feat)))
dat <- merge(dat,
             feat, 
             by.x= "L",
             by.y= "ID")
dat <- merge(dat,
             feat, 
             by.x= "R",
             by.y= "ID",
             suffixes= c("_L", "_R"))
# formula
.f <- paste(c(paste0("`", names(feat)[-1], "_L`"), paste0("`", names(feat)[-1], "_R`")), collapse= "+")

# Define train and test sets
set.seed(1)
dat[dat[, .(set= sample(5)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(5)), R], setR:= i.set, on= "R"]
dat[, set:= paste0("test", .GRP), .(setL, setR)]
# Train linear model for each train set and compute predicted values
dat[, c("predicted", "R2"):= {
  print(set)
  cL <- L
  cR <- R
  train <- dat[!(L %in% cL) & !(R %in% cR)]
  model <- lm(formula = as.formula(paste0("residuals~median_L*median_R+", .f)),
              data= train)
  pred <- predict(model, newdata = .SD)
  rsq <- vl_model_eval(observed = residuals, 
                       predicted = pred)$Rsquare
  .(pred, rsq)
}, set]

dat[, {
  smoothScatter(residuals, predicted)
  legend("topleft",
         paste0("PCC= ", round(cor.test(residuals, predicted)$estimate, 2)),
         bty= "n")
}]

saveRDS(dat[, .(L, R, predicted, R2)], 
        "Rdata/CV_linear_model_residuals_vllib002.rds")


