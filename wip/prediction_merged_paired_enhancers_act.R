setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
options(datatable.print.topn=1)
require(caret)

# #--------------------------------------------------------------#
# # 1- Format data
# #--------------------------------------------------------------#
# dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
# dat <- dat[!is.na(median_L) & !is.na(median_R) & !is.na(log2FC_add) & !is.na(log2FoldChange) & !is.na(diff) & enh_L!=enh_R]
# 
# feat <- readRDS("Rdata/library/lib_features.rds")
# cols <- c("ID", grep("^motif__", colnames(feat), value = T), "group")
# feat <- feat[, ..cols]
# 
# obj <- copy(dat)
# obj <- obj[, .(enh_L, enh_R, 
#                log2FoldChange= factor(ifelse(log2FoldChange>2, "yes", "no")), 
#                super= factor(ifelse(diff>1, "yes", "no")), 
#                sub= factor(ifelse(diff< -1, "yes", "no")))]
# cols <- c("ID", grep("^motif", colnames(feat), value = T))
# obj_L <- merge(obj, feat[, ..cols], by.x= "enh_L", by.y= "ID")
# obj_R <- merge(obj, feat[, ..cols], by.x= "enh_R", by.y= "ID")
# obj <- rbind()
# obj <- rbind(melt(obj_L, id.vars = c("enh_L", "enh_R", "log2FoldChange", "super", "sub")),
#              melt(obj_R, id.vars = c("enh_L", "enh_R", "log2FoldChange", "super", "sub")))
# obj <- obj[, .(value= sum(value)), .(enh_L, enh_R, log2FoldChange, super, sub, variable)]
# obj <- dcast(obj, enh_L+enh_R+log2FoldChange+super+sub~variable, value.var = "value")
# 
# #--------------------------------------------------------------#
# # 2- Test models
# #--------------------------------------------------------------#
# cmd <- "module load r/3.6.2-foss-2018b; Rscript /groups/stark/vloubiere/pipelines/classification_prediction_models.R"
# out_fold <- "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/motifs/models/merged_enhancers/"
# for(class_var in c("log2FoldChange", "super", "sub"))
# {
#   for(model in c("logistic", "rpart", "rf", "nb", "nb_kernel", "svmLinear2", "svmRadial"))
#   {
#     for(sample in c("class", "dummy"))
#     {
#       output <- paste0(out_fold, class_var, "_", sample, "_prediction_", model, "_model.rds")
#       cols <- c(grep("^motif", colnames(obj), value = T), class_var)
#       tmp <- tempfile(fileext = ".rds")
#       current <- obj[, ..cols]
#       if(sample=="dummy")
#       {
#         current[[class_var]] <- sample(current[[class_var]], nrow(current))
#       }
#       saveRDS(current, tmp)
#       bsub(cmd= paste(cmd, tmp, class_var, model, output), m = 20, o = out_fold, e= out_fold, t= "2-00:00:00")
#       print(paste(class_var, model, sample, "  --->  SUBMITTED!"))
#     }
#   }
# }

dat <- data.table(file= list.files("Rdata/motifs/models/merged_enhancers/", ".rds", full.names = T))
dat[, c("class_var", "sample"):= tstrsplit(basename(file), "_|[.]", keep= c(1,2)), (dat)]
dat <- dat[, readRDS(file)$diagnostic, .(file, class_var, sample)]
res <- dcast(dat, class_var+sample~model, value.var = "PR.AUC")

par(mfrow= c(1,3), las=1, mar= c(5,5,10,5))
res[, 
    {
      current <- copy(.SD)
      setkeyv(current, "sample")
      barplot(as.matrix(current[c("dummy", "class"), logistic:rpart]), beside= T, ylab= "Precision Recall AUC",
              main= class_var[1], col= c("darkgrey", "lightgrey"), ylim= c(0, 1.2))
      legend("topleft", fill= c("darkgrey", "lightgrey"), legend = c("random", "trained"))
      print("")
    }, class_var]





