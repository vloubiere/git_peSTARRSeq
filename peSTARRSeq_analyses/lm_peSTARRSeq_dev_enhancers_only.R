setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
options(datatable.print.topn= 1)
require(data.table)
require(pheatmap)
require(parallel)

#------------------------------------------------------------------#
# 1- Format data
#------------------------------------------------------------------#
# Import data
dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[!is.na(median_L) & !is.na(median_R) & !is.na(log2FoldChange) & !is.na(diff)]
rep <- readRDS("Rdata/processed_peSTARRSeq_data/DSCPr1_inputr2_DSCPr2_inputr3.rds")
dat[rep, rep1:= i.log2FoldChange, on= c("enh_L", "enh_R")]
rep <- readRDS("Rdata/processed_peSTARRSeq_data/DSCPr3_inputr4_DSCPr4_inputr5.rds")
dat[rep, rep2:= i.log2FoldChange, on= c("enh_L", "enh_R")]
# Add motif counts
feat <- readRDS("Rdata/library/lib_features.rds")
cols <- c("ID", grep("^motif", colnames(feat), value= T))
feat <- feat[, ..cols]
setkeyv(feat, "ID")
dat <- rbind(cbind(dat, feat[dat$enh_L, !"ID"]),
             cbind(dat, feat[dat$enh_R, !"ID"]))
dat <- dat[, lapply(.SD, sum), ID:rep2]
dat[, group_L:= tstrsplit(enh_L, "_", keep=1)]
dat[, group_R:= tstrsplit(enh_R, "_", keep=1)]

# Select only developmental enhancer pairs
dat <- dat[grepl("^dev", enh_L) & grepl("^dev", enh_R)]

#------------------------------------------------------------------#
#2- Modeling
#------------------------------------------------------------------#
# motif_pairwise_cmb <- CJ(grep("^motif", colnames(dat), value = T),
#                          grep("^motif", colnames(dat), value = T))
# motif_pairwise_cmb[, cmb:= paste(c(V1, V2), collapse = "*"), (motif_pairwise_cmb)]
# pred <- data.table(form= c("median_L",
#                            "median_R",
#                            "median_L+median_R",
#                            "median_L*median_R*log2FC_add",
#                     paste0("median_L*median_R*log2FC_add+", paste(motif_pairwise_cmb$cmb, collapse = "+"))),
#                    name= c("Left individual activity",
#                            "Right individual activity",
#                            "Left+Right individual activities",
#                            "Left*Right*Exp.Additive",
#                            "Left*Right*Exp.Additive + pairwise motif combinations w/ interaction"))
# pred[, form:= paste0("log2FoldChange~", form), form]
# models <- list()
# for(i in seq(nrow(pred)))
# {
#   models[[pred[i, name]]] <- lm(as.formula(pred[i, form]), dat)
#   print(i)
# }
# saveRDS(models, "Rdata/modeling_peSTARRSeq/linear_models_dev_enhancers_only.rds")

#------------------------------------------------------------------#
# 3- Plots
#------------------------------------------------------------------#
if(!exists("mod"))
{
  mod <- readRDS("Rdata/modeling_peSTARRSeq/linear_models_dev_enhancers_only.rds")
}


pdf("pdf/peSTARRSeq/log2FC_lm_w_wo_interctions_dev_enhancers_only.pdf", 9, 10)
par(mfrow= c(3,3), las= 1)
for(i in seq(mod))
{
  smoothScatter(predict(mod[[i]]), dat$log2FoldChange, xlab= "Predicted value", ylab= "log2FoldChange")
  mtext(names(mod)[i], side= 3, line = 1, cex= 0.6)
  abline(0, 1)
  rsq <- paste("R²=", round(summary(mod[[i]])$r.squared, 2))
  PCC <- paste("PCC=", round(cor.test(predict(mod[[i]]), dat$log2FoldChange)$estimate, 2))
  legend("topleft", legend= c(rsq, PCC), bty= "n")
  #my_fig_label(c("A", "B", "C", "D", "E", "F", "G", "H")[i], cex=2)
}

smoothScatter(dat$rep1, dat$rep2, xlab= "act. replicates 1-2", ylab= "act. replicates 3-4")
mtext("Replicates", side= 3, line = 1, cex= 0.6)
abline(0, 1)
PCC <- paste("PCC=", round(cor.test(dat$rep1, dat$rep2)$estimate, 2))
legend("topleft", legend= PCC, bty= "n")
# my_fig_label("I", cex=2)
dev.off()



