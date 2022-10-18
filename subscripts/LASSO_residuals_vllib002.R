setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- unique(dat[actClassL!= "inactive", .(L, indL)])[indL<=quantile(indL, 0.9), L]
QR <- unique(dat[actClassR!= "inactive", .(R, indR)])[indR<=quantile(indR, 0.9), R]
dat <- dat[L %in% QL & R %in% QR]

dat <- rbind(unique(dat[, .(side= "5'", ID= L, ind= indL, meanResiduals= meanResidualsL)]),
             unique(dat[, .(side= "3'", ID= R, ind= indR, meanResiduals= meanResidualsR)]))
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
sel <- vl_Dmel_motifs_DB_full[!is.na(FBgn), motif_ID]
counts <- vl_motif_counts(lib[dat, enh_sequence, on= "ID_vl==ID"], sel= sel)
dat <- cbind(dat[, .(side, ID, meanResiduals, ind)], counts)

# LASSO regression
model <- dat[, {
  x <- as.matrix(.SD[, !"ID"][,-1])
  y <- as.matrix(.SD[, !"ID"][,1])
  # perform k-fold cross-validation to find optimal lambda value
  set.seed(1)
  cv_model <- cv.glmnet(as.matrix(x), y, alpha = 1)
  # find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  # Model
  set.seed(1)
  best_model <- glmnet(x= as.matrix(x), 
                       y= y, 
                       alpha = 1, 
                       lambda = best_lambda)
  # Assess performance
  predicted <- predict(best_model, newx = as.matrix(x))
  rsq <- vl_model_eval(y, predicted)$Rsquare
  # Return
  best_model$predict <- .SD[, .(ID, meanResiduals, predicted= predicted)]
  best_model$coeffs <- as.data.table(as.matrix(coef(best_model)), 
                                     keep.rownames = var)
  best_model$rsq <- rsq
  .(model= .(best_model))
}, side]

modL <- model[side=="5'", model[[1]]]
coeffs <- modL$coeffs[s0!=0][order(s0)]
coeffs[vl_Dmel_motifs_DB_full, name:= i.motif_cluster, on= "rn==motif_ID"]
coeffs[rn=="(Intercept)", name:= "Intercept"]
coeffs <- coeffs[, .SD[which.max(abs(s0))], name]

########################
# CV sub models
########################
CV_rsq <- mclapply(1:25, function(i)
{
  set.seed(i)
  sample <- sample(length(y), round(length(y)*0.9))
  cv_model <- cv.glmnet(as.matrix(x)[sample,], y[sample], alpha = 1)
  best_lambda <- cv_model$lambda.min
  best_model <- glmnet(x= as.matrix(x)[sample,], 
                       y= y[sample], 
                       alpha = 1, 
                       lambda = best_lambda)
  predicted <- predict(best_model, newx = as.matrix(x)[-sample,])
  print(i)
  # Return
  vl_model_eval(y[-sample], predicted)$Rsquare
}, mc.preschedule = F, mc.cores = getDTthreads())