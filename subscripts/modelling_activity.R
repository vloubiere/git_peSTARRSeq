# Smoothscatter models ####
if(!file.exists("Rdata/linear_models_prediction.rds"))
{
  feat <- readRDS("Rdata/master_lib_features.rds")
  cols <- c("ID", grep("^motif", colnames(feat), value= T))
  feat <- feat[, ..cols]
  tab <- readRDS("Rdata/master_results_peSTARRSeq.rds")
  tab <- tab[!(spike_in)]
  tab <- feat[tab, , on= "ID==L"]
  tab <- merge(tab, feat, by.x= "R", by.y= "ID", all.x= T, suffixes= c("___L", "___R"))
  cols <- unique(gsub("___L|___R", "", grep("^motif", colnames(tab), value = T)))
  tab[, (paste0(cols, "___sum")):= lapply(cols, function(x) rowSums(tab[, grep(x, colnames(tab)), with= F], na.rm= T))]
  mots1mots2 <- paste0(c(grep("___L", colnames(tab), value = T), grep("___R", colnames(tab), value = T)), collapse= "+")
  motssum <- paste0(grep("___sum", colnames(tab), value = T), collapse= "+")
  motspw <- CJ(grep("___L", colnames(tab), value = T), grep("___R", colnames(tab), value = T))[, paste0(V1, "*", V2)]
  motspw <- paste0(motspw, collapse= "+")
  models <- data.table(lib= unique(tab$lib))
  models <- models[, .(title= c("lm_L", 
                                "lm_R", 
                                "lm_add", 
                                "lm_L+R", 
                                "lm_L*R", 
                                "lm_L*R+add", 
                                "lm_L*R*add", 
                                "lm_mots1+mot2", 
                                "lm_motsums", 
                                "lm_motspw", 
                                "lm_L+R+mots1+mot2", 
                                "lm_L+R+motsums", 
                                "lm_L+R+motspw")), lib]
  models[, file:= paste0("db/models/", lib, "_", title, ".rds"), models]
  if(any(!file.exists(models$file)))
  {
    models[, formula:= c("log2FoldChange~median_L",
                         "log2FoldChange~median_R",
                         "log2FoldChange~add",
                         "log2FoldChange~median_L+median_R",
                         "log2FoldChange~median_L*median_R",
                         "log2FoldChange~median_L*median_R+add",
                         "log2FoldChange~median_L*median_R*add",
                         paste0("log2FoldChange~", mots1mots2),
                         paste0("log2FoldChange~", motssum),
                         paste0("log2FoldChange~", motspw),
                         paste0("log2FoldChange~median_L+median_R+", mots1mots2),
                         paste0("log2FoldChange~median_L+median_R+", motssum),
                         paste0("log2FoldChange~median_L+median_R+", motspw)), lib]
    setkeyv(tab, "lib")
    models[, {if(!file.exists(file)){saveRDS(lm(as.formula(formula), tab[lib]), file)}}, models]
    models$formula <- NULL
  }
  models[, rsq:= summary(readRDS(file))$r.squared, file]
  models[, obs:= .(list(readRDS(file)$model$log2FoldChange)), file]
  models[, exp:= .(list(predict(readRDS(file)))), file]
  models[, PCC:= cor.test(unlist(obs), unlist(exp))$estimate, .(lib, title)]
  saveRDS(models, "Rdata/linear_models_prediction.rds")
}else{
  models <- readRDS("Rdata/linear_models_prediction.rds")
}

pdf("pdf/linear_models_peSTARRSeq.pdf", width = 52, height = 4)
par(mfrow= c(1,13), las= 1)
models[, {
  smoothScatter(unlist(exp), unlist(obs), main= paste(lib, title))
  legend("topleft", c(paste0("R2= ", round(rsq, 2)), paste0("PCC= ", round(PCC, 2))), bty= "n")
  print("")
}, .(lib, title, PCC, rsq)]
dev.off()
####

