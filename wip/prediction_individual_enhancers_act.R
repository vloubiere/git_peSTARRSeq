setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
options(datatable.print.topn=1)
require(kohonen)
require(caret)

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
obj[, median_L:=  factor(ifelse(median_L>2, "yes", "no"))]
obj[, median_R:=  factor(ifelse(median_R>2, "yes", "no"))]

#--------------------------------------------------------------#
# 2- Test models
#--------------------------------------------------------------#
### Left
if(!file.exists("Rdata/motifs/prediction_ind_activities.R"))
{
  # Dummy pred
  dummy <- copy(obj)
  dummy[, median_L:= sample(median_L, nrow(obj))]
  dummy[, median_R:= sample(median_R, nrow(obj))]
  cols <- c(grep("^motif", colnames(obj), value = T), "median_L")
  medL_dummy <- my_mult_models(dummy[, ..cols], class_var = "median_L")
  cols <- c(grep("^motif", colnames(obj), value = T), "median_R")
  medR_dummy <- my_mult_models(dummy[, ..cols], class_var = "median_R")
  
  # pred
  cols <- c(grep("^motif", colnames(obj), value = T), "median_L")
  medL_pred <- my_mult_models(obj[, ..cols], class_var = "median_L")
  cols <- c(grep("^motif", colnames(obj), value = T), "median_R")
  medR_pred <- my_mult_models(obj[, ..cols], class_var = "median_R")
  
  all.res <- list(medL_dummy= medL_dummy, medR_dummy= medR_dummy, 
                  medL_pred= medL_pred, medR_pred= medR_pred)
  saveRDS(all.res, "Rdata/motifs/prediction_ind_activities.R")
}
all.res <- readRDS("Rdata/motifs/prediction_ind_activities.R")

res <- list(ctl_L= rbindlist(all.res$medL_dummy$diag)$PR.AUC,
            act_L= rbindlist(all.res$medL_pred$diag)$PR.AUC,
            ctl_R= rbindlist(all.res$medR_dummy$diag)$PR.AUC,
            act_R= rbindlist(all.res$medR_pred$diag)$PR.AUC)

pdf("pdf/peSTARRSeq/prediction_individual_enhancer_activities.pdf", width = 10)
par(mfrow= c(1, 2))
pl <- rbind(res$ctl_L, res$act_L)
colnames(pl) <- names(all.res$medL_dummy$trained.models)
pl <- pl[, order(pl[2,])]
bar <- barplot(pl, beside = T, ylim= c(0, 1), las= 2, names.arg = rep("", ncol(pl)), main= "Left position",
               ylab= "Precision recall AUC")
legend("topleft", legend = c("Random", "Precision recall AUC"), bty= "n", fill= c("black", "lightgrey"))
text(colMeans(bar), -0.03, colnames(pl), srt= 45, xpd= T, pos= 2, offset= -0.25)

pl <- rbind(res$ctl_R, res$act_R)
colnames(pl) <- names(all.res$medL_dummy$trained.models)
pl <- pl[, order(pl[2,])]
bar <- barplot(pl, beside = T, ylim= c(0, 1), las= 2, names.arg = rep("", ncol(pl)), main= "Right position",
               ylab= "Precision recall AUC")
legend("topleft", legend = c("Random", "Precision recall AUC"), bty= "n", fill= c("black", "lightgrey"))
text(colMeans(bar), -0.03, colnames(pl), srt= 45, xpd= T, pos= 2, offset= -0.25)
dev.off()
