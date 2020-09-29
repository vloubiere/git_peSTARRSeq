setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
options(datatable.print.topn=1)
require(kohonen)
require(randomForest)

#--------------------------------------------------------------#
# 1- Format data
#--------------------------------------------------------------#
dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[!is.na(median_L) & !is.na(median_R) & !is.na(log2FC_add) & !is.na(log2FoldChange) & !is.na(diff) & enh_L!=enh_R]
feat <- readRDS("Rdata/library/lib_features.rds")

enh_L <- unique(dat[, .(ID= enh_L, median_L)])
enh_R <- unique(dat[, .(ID= enh_R, median_R)])
cols <- c("ID", grep("^motif__", colnames(feat), value = T), "group")
obj <- feat[ID %in% dat$enh_L | ID %in% dat$enh_R, ..cols]
obj <- obj[enh_L, , on= "ID"]
obj <- obj[enh_R, , on= "ID"]
obj <- na.omit(obj)

#--------------------------------------------------------------#
# 2- Define train and test sets
#--------------------------------------------------------------#
dat <- copy(obj)
dat <- dat[, act_L:= cut(median_L, c(-Inf, 1, 4, 7, Inf), labels = c("inactive", "weak", "medium", "strong"))]
dat <- dat[, act_R:= cut(median_R, c(-Inf, 1, 4, 7, Inf), labels = c("inactive", "weak", "medium", "strong"))]
set.seed(1)
sel <- sample(nrow(obj), nrow(obj)*0.7)
train <- dat[sel]
test <- dat[-sel]

res <- list(train_set= train, test_set= test, 
            models= NULL, pred= NULL, accuracy= NULL)

#--------------------------------------------------------------#
# 3- SOM
#--------------------------------------------------------------#
som.train <- list(motif= as.matrix(train[, ID:motif__predrem__nrMotif1321], 1),
                  group= factor(train[, group]),
                  act_L= factor(train[, act_L]),
                  act_R= factor(train[, act_R]))
som.test <- list(motif= as.matrix(test[, ID:motif__predrem__nrMotif1321], 1),
                 group= factor(test[, group]),
                 act_L= factor(test[, act_L]),
                 act_R= factor(test[, act_R]))

# Training
mygrid <- somgrid(xdim= 4, ydim= 4, topo = 'hexagonal', toroidal = T)
set.seed(1234)
som <- supersom(som.train, grid = mygrid, rlen = 500, maxNA.fraction = 0.3)
res$models[["som"]] <- som

# Prediction
pred <- predict(som, newdata = som.test, whatmap = 1)
accuracy <- lapply(c("group", "act_L", "act_R"), function(x)
{
  res <- as.data.table(table(som.test[[x]], pred$predictions[[x]]))
  res[, cdition:= x]
  return(res)
})
accuracy <- rbindlist(accuracy)
accuracy <- accuracy[, sum(.SD[V1==V2, N])/sum(.SD[, N])*100, cdition]

res$pred[["som"]] <- pred
res$accuracy[["som"]] <- accuracy

#--------------------------------------------------------------#
# 4- Random forest
#--------------------------------------------------------------#
for(var in c("group", "act_L", "act_R"))
{
  current <- train[, !c("ID", "group", "median_L", "median_R", "act_L", "act_R")]
  current <- cbind(current, train[, lapply(.SD, factor), .SDcols= var])
  .f <- as.formula(paste0(var, "~."))
  rf <- randomForest(formula= .f, data= current, ntree=100, mtry= round(sqrt(ncol(current)-1)), importance= TRUE)
  
  cols <- setdiff(colnames(current), var)
  pred <- predict(rf, newdata= test[, ..cols])
  accuracy <- length(which(pred==test[[var]]))/nrow(test)*100
  
  res$models[["rf"]][[var]] <- rf
  res$pred[["rf"]][[var]] <- pred
  res$accuracy[["rf"]][[var]] <- accuracy
}

#--------------------------------------------------------------#
# 5- Linear models
#--------------------------------------------------------------#
for(var in c("median_L", "median_R"))
{
  current <- train[, !c("ID", "group", "median_L", "median_R", "act_L", "act_R")]
  current <- cbind(current, train[, ..var])
  .f <- as.formula(paste0(var, "~."))
  .lm <- lm(formula= .f, data= current)
  res$models[["lm"]][[var]] <- .lm
  
  cols <- setdiff(colnames(current), var)
  pred <- predict(.lm, newdata= test[, ..cols])
  accuracy <- summary(lm(pred~test[[var]]))$r.squared
  
  res$models[["lm"]][[var]] <- rf
  res$pred[["lm"]][[var]] <- pred
  res$accuracy[["lm"]][[var]] <- accuracy
}

saveRDS(res, "Rdata/motifs/prediction_individual_act.rds")





