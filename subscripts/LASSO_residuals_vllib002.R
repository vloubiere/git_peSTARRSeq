setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh."]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
lib[, predicted:= predict(model, newdata = .SD)]
lib[, residuals:= log2FoldChange-predicted]
feat <- fread("Rdata/final_300bp_enhancer_features.txt")

# Train models to predict left and right median residuals
if(!file.exists("Rdata/LASSO_residuals_vllib002.rds"))
{
  models <- list()
  for(pos in c("L", "R"))
  {
    # Retrieve data
    .c <- lib[, .(median_residuals= median(residuals)), c(pos, paste0("median_", pos))]
    setnames(.c, c("ID", "median", "median_residuals"))
    .c <- cbind(.c, feat[.c, on= "ID"])
    
    # Split predictors and dependent var
    x <- cbind(.c[, .(median)],
               # .c[, .(ATAC, H3K27Ac, H3K4me1, H3K4me3, H3K27me3, GAF, SUHW)], K4me1 and 3 are informative, but not top predictors
               .c[, names(.c) %in% vl_Dmel_motifs_DB_full$motif, with= F])
    y <- .c$median_residuals
    
    ########################
    # Full model
    ########################
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
    # Return
    rsq <- vl_model_eval(y, predicted)$Rsquare
    
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
    
    models[[pos]] <- list(model= best_model,
                          predict= .c[, .(ID, median_residuals, predicted= predicted)],
                          coeffs= as.data.table(as.matrix(coef(best_model)), 
                                                keep.rownames = var),
                          rsq= rsq,
                          CV_rsq= unlist(CV_rsq))
  }
  
  saveRDS(models,
          "Rdata/LASSO_residuals_vllib002.rds")
}else
  models <- readRDS("Rdata/LASSO_residuals_vllib002.rds")

######################################################
# Select variables of interest
######################################################
top <- rbind(models$L$coeffs, models$R$coeffs)
top <- top[vl_Dmel_motifs_DB_full, motif_cluster:= motif_cluster, on= "rn==motif"]
top <- na.omit(top)[order(abs(s0), decreasing = T)][s0!=0]
top <- top[, .SD[1], motif_cluster]
top <- top[, lib[, .(L, R, median_L, median_R, residuals)], rn]
motifs <- c("DRE/1__2", "lola/10__6", "AP1/1__42", "Ebox/CACGTG/1__1", "D19A__1",
            "TBX/1__13", "Ebox/CATATG/twi__2", "GCM/1__3", "Mad/1__1", "Ebox/CAGCTG/3__14",
            "CG11152/7", "NFAT/1__8", "HSF/4__2", "Ebox/CACGTG/3__42", "KLF/SP/1__19", "tgo", 
            "CG32830/br/vvl__56")
top <- top[motifs, on= "rn"]
feats <- melt(feat, 
              id.vars = "ID", 
              measure.vars = unique(top$rn))
top[feats, cut_L:= i.value, on= c("L==ID", "rn==variable")]
top[feats, cut_R:= i.value, on= c("R==ID", "rn==variable")]
cols <- c("cut_L", "cut_R")
top[, (cols):= lapply(.SD, function(x) {
  .q <- quantile(x, 0.995)
  ifelse(x>.q, .q, x)
}), rn, .SDcols= cols]
top[, cut:= paste0(cut_L, "_", cut_R)]
coeffs <- top[, {
  .c <- summary(lm(residuals~median_L*median_R+cut))$coefficients
  .c <- as.data.table(.c, keep.rownames = T)[grepl("^cut", rn)]
  .c[, c("L", "R") := lapply(tstrsplit(gsub("^cut", "", rn), "_"), as.integer)]
  .(.(as.matrix(dcast(.c, L~R, value.var = "Estimate"), 1)))
}, rn]

######################################################
# PLOT
######################################################
pdf("pdf/draft/lasso_predictors_interaction_residuals.pdf", 4, 4)
par(mar= c(4, 4, 4, 4),
    tcl= -0.2,
    las= 1)
lib[, strength_L:= cut(median_L, c(-Inf, 2,4,6, Inf), c("<2", "2-4", "4-6", ">6"))]
lib[, strength_R:= cut(median_R, c(-Inf, 2,4,6, Inf), c("<2", "2-4", "4-6", ">6"))]
mat <- dcast(lib, strength_L~strength_R, value.var = "residuals", fun.aggregate = mean)
mat <- mat[nrow(mat):1,]
vl_heatmap(as.matrix(mat, 1), 
           cluster_rows = F, 
           cluster_cols = F,
           breaks= c(-1,0,1),
           legend_title = "Obs./Exp. (log2)",
           main= "Activity")
coeffs[, {
  .c <- V1[[1]]
  vl_heatmap(.c[nrow(.c):1,],
             cluster_rows = F,
             cluster_cols = F,
             breaks= c(-1,0,1),
             main= rn,
             show_legend= F)
  vl_heatkey(breaks = c(-1,0,1), 
             left= par("usr")[2]+strwidth("M"),
             col= c("cornflowerblue", "white", "red"), 
             main = "lm coeff.")
  print("")
}, rn]
dev.off()

