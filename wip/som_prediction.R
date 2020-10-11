setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
options(datatable.print.topn=1)
require(kohonen)

#--------------------------------------------------------------#
# 1- Train som using motifs on a subset of the data
#--------------------------------------------------------------#
dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[!is.na(median_L) & !is.na(median_R) & !is.na(log2FC_add) & !is.na(log2FoldChange) & !is.na(diff) & enh_L!=enh_R]
feat <- readRDS("Rdata/library/lib_features.rds")

enh_L <- unique(dat[, .(ID= enh_L, act_L= cut(median_L, c(-Inf, 1, 4, 7, Inf), labels = c("inactive", "weak", "medium", "strong")))])
enh_R <- unique(dat[, .(ID= enh_R, act_R= cut(median_R, c(-Inf, 1, 4, 7, Inf), labels = c("inactive", "weak", "medium", "strong")))])
enh_L <- unique(dat[, .(ID= enh_L, act_L= median_L)])
enh_R <- unique(dat[, .(ID= enh_R, act_R= median_R)])
cols <- c("ID", grep("^motif__", colnames(feat), value = T), "group")
obj <- feat[ID %in% dat$enh_L | ID %in% dat$enh_R, ..cols]
obj <- obj[enh_L, , on= "ID"]
obj <- obj[enh_R, , on= "ID"]

sel <- sample(nrow(obj), nrow(obj)*0.7)
train <- obj[sel]
train <- list(motif= as.matrix(train[, ID:motif__predrem__nrMotif1321], 1),
              group= factor(train[, group]),
              act_L= factor(train[, act_L]),
              act_R= factor(train[, act_R]))

mygrid <- somgrid(xdim= 4, ydim= 4, topo = 'hexagonal', toroidal = T)
set.seed(1234)
som <- supersom(train, grid = mygrid, rlen = 500, maxNA.fraction = 0.3)

#--------------------------------------------------------------#
# 2- Check predictive power of the som on enhancer group and left right ind act
#--------------------------------------------------------------#
test <- obj[-sel]
test <- list(motif= as.matrix(test[, ID:motif__predrem__nrMotif1321], 1),
             group= factor(test[, group]),
             act_L= factor(test[, act_L]),
             act_R= factor(test[, act_R]))

pred <- predict(som, newdata = test, whatmap = 1)
accuracy <- lapply(c("group", "act_L", "act_R"), function(x)
{
  res <- as.data.table(table(test[[x]], pred$predictions[[x]]))
  res[, cdition:= x]
  return(res)
})
accuracy <- rbindlist(accuracy)
accuracy <- accuracy[, sum(.SD[V1==V2, N])/sum(.SD[, N])*100, cdition]

pdf("pdf/peSTARRSeq/som_individual_act_prediction_motifs.pdf", width = 9)
par(mfrow= c(2,3))
plot(som, "counts", shape= "straight")
plot(som, "codes", shape= "straight")
par(mar= c(5,5,5,5), las= 1)
barplot(accuracy$V1, names.arg = accuracy$cdition, main= "Group", ylab= "percentage accuracy")
dev.off()
