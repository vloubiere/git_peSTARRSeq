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
cols <- c("ID", grep("^motif__", colnames(feat), value = T), "group")
feat <- feat[, ..cols]

dat[, act:= cut(log2FoldChange, c(-Inf, 2, 5, 3, Inf), labels = c("inactive", "weak", "medium", "strong"))]
dat <- merge(dat[, .(ID, enh_L, enh_R, log2FoldChange, diff)], feat, by.x= "enh_L", by.y= "ID")
dat <- merge(dat, feat, by.x= "enh_R", by.y= "ID", suffixes= c("_L", "_R"))
.m <- melt(dat, id.vars = c("ID", "enh_L", "enh_R", "group_L", "group_R", "log2FoldChange", "diff"), 
           measure.vars = patterns("^motif"))
.m[, variable:= gsub("_L|_R", "", variable)]
dmat <- .m[, .(value= sum(value)), setdiff(colnames(.m), "value")]
dmat <- dcast(dmat, ID+enh_L+enh_R+group_L+group_R+log2FoldChange+diff~variable, value.var = "value")
#--------------------------------------------------------------#
# 2- Define train and test sets
#--------------------------------------------------------------#
set.seed(1)
sel <- sample(nrow(dmat), nrow(dmat)*0.7)
train <- dmat[sel]
test <- dmat[-sel]

res <- list(train_set= train, test_set= test, 
            models= NULL, pred= NULL, accuracy= NULL)

#--------------------------------------------------------------#
# 3- SOM
#--------------------------------------------------------------#
cols <- c("ID", grep("^motif", colnames(train), value= T))
som.train <- list(motif= as.matrix(train[, ..cols], 1),
                  group_L= factor(train[, group_L]),
                  group_R= factor(train[, group_R]),
                  act= factor(cut(train[, log2FoldChange], c(-Inf, 1, 4, 7, Inf), c("inactive", "weak", "medium", "strong"))),
                  diff= factor(cut(train[, diff], c(-Inf, -1.5, 1.5, Inf), c("sub-additive", "additive", "super-additive"))))
som.test <- list(motif= as.matrix(test[, ..cols], 1),
                 group_L= factor(test[, group_L]),
                 group_R= factor(test[, group_R]),
                 act= factor(cut(test[, log2FoldChange], c(-Inf, 1, 4, 7, Inf), c("inactive", "weak", "medium", "strong"))),
                 diff= factor(cut(test[, diff], c(-Inf, -1.5, 1.5, Inf), c("sub-additive", "additive", "super-additive"))))
# Training
mygrid <- somgrid(xdim= 18, ydim= 18, topo = 'hexagonal', toroidal = T)
set.seed(1234)
som <- supersom(som.train, grid = mygrid, rlen = 250, maxNA.fraction = 0.3)
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





