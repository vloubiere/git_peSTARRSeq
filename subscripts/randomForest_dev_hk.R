setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(randomForest)


#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/final_results_table.rds")
# Restrict to "smaller" and cleaner twist 12
dat <- dat[library=="T12"]
setkeyv(dat, c("CP", "spacer", "group_L", "group_R"))
# Select dependent & predictor variables
obj <- rbindlist(list(dev= dat[.("DSCP", "SCR1_300")],
                      hk= dat[.("RpS12", "SCR1_300")]), 
                 idcol = "group")
# Add motifs
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
sel_motifs <- unique(feat$motif_padj[padj<1e-5 & log2OR>0, motif])
sel_motifs <- as.character(sel_motifs)
obj <- feat$add_feature(DT= obj, 
                        feature = feat$top_motifs[, c("ID", sel_motifs), with= F])
# Add chromatin features and replace NAs with 0
obj <- feat$add_feature(DT= obj, 
                        feature = feat$chromatin_features)
chrom <- c(paste0(names(feat$chromatin_features)[-1], "_L"),
           paste0(names(feat$chromatin_features)[-1], "_R"))
obj[, (chrom):= lapply(.SD, function(x) {x[is.na(x)] <- 0; return(x)}), .SDcols= chrom]

#-----------------------------------------------#
# MODELLING
#-----------------------------------------------#
# Training variables
vars <- c("median_L", "median_R", 
          paste0(sel_motifs, "_L"), 
          paste0(sel_motifs, "_R"),
          chrom)

#-----------------------------------------------#
# randomForest
#-----------------------------------------------#
models <- obj[, {
  # Define train and test subsets
  seed <- 1
  set.seed(seed)
  selL <- sample(unique(L), 0.85*length(unique(L)))
  set.seed(seed)
  selR <- sample(unique(L), 0.85*length(unique(L)))
  train <- L %in% selL & R %in% selR
  
  data <- as.matrix(.SD)
  X_train <- data[train,]
  Y_train <- log2FoldChange[train]
  X_test <- data[!train,]
  Y_test <- log2FoldChange[!train]
  
  # Model
  model <- randomForest(x = X_train, 
                        y= Y_train,
                        xtest= X_test,
                        ytest= Y_test,
                        importance= T)
  
  # RETURN
  .(model= .(model))
}, group, .SDcols= vars]

saveRDS(models,
        "test/randomForest.rds")