setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
dir.create("db/linear_models", 
           showWarnings = F)
dir.create("db/linear_models/global_models", 
           showWarnings = F)
dir.create("db/linear_models/motif_pairs", 
           showWarnings = F)

# Import
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")
dat <- DT[L!=R]
dat <- dat[!is.na(additive)]
dat <- dat[(cdition== "vllib002" & group_L== "dev" & group_R== "dev" & median_L>1 & median_R>1) |
           (cdition== "vllib016" & group_L== "hk" & group_R== "hk")]
# Top motifs
comb <- data.table(L_mot= grep("^motif__.*_L", names(dat), value = T))
comb[, R_mot:= gsub("_L$", "_R", L_mot)]
comb[, colName:= gsub("_L$", "_merge", L_mot)]
comb[, {
  dat[[colName]] <<- dat[[.BY[[1]]]]+dat[[.BY[[2]]]]
}, (comb)]
# Som motif groups
comb <- data.table(L_mot= grep("^motifSom__.*_L", names(dat), value = T))
comb[, R_mot:= gsub("_L$", "_R", L_mot)]
comb[, colName:= gsub("_L$", "_merge", L_mot)]
comb[, {
  dat[[colName]] <<- dat[[.BY[[1]]]]+dat[[.BY[[2]]]]
}, (comb)]
setkeyv(dat, "active_plot_group")

#-------------------------#
# Build models formulas
#-------------------------#
# Make a string for all motifs
mot <- CJ(grep("^motif__.*_L$", names(dat), value = T), grep("^motif__.*_R$", names(dat), value = T), unique= T)
motSom <- CJ(grep("^motifSom__.*_L$", names(dat), value = T), grep("^motifSom__.*_R$", names(dat), value = T), unique= T)
motMerge <- grep("^motif__.*_merge$", names(dat), value = T)
motMergeSom <- grep("^motifSom__.*_merge$", names(dat), value = T)
mod <- rbind(data.table(form= "log2FoldChange~median_L*median_R", # log2FC / Activity only w/ interaction
                        save= "predict_log2FoldChange__actL*R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~median_L*median_R+", paste0(paste0(mot[, paste0(V1, "*", V2)], collapse= "+"))), # log2FC / Activity and motifs w/ interaction
                        save= "predict_log2FoldChange__actL*R+motL*R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~median_L*median_R+", paste0(paste0(motSom[, paste0(V1, "*", V2)], collapse= "+"))), # log2FC / Activity and motifsSOM w/ interaction
                        save= "predict_log2FoldChange__actL*R+motSomL*R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~", paste0(paste0(mot[, paste0(V1, "*", V2)], collapse= "+"))),  # log2FC / Motifs only w/ interaction
                        save= "predict_log2FoldChange__motL*R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~", paste0(paste0(motSom[, paste0(V1, "*", V2)], collapse= "+"))),  # log2FC / MotifsSOM only w/ interaction
                        save= "predict_log2FoldChange__motSomL*R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~", mot[, paste0(V1, "+", V2)]), # log2FC / Motif pair w/0 interaction
                        save= paste0("predict_log2FoldChange__", mot[, paste0(V1, "+", V2)]),
                        folder= "db/linear_models/motif_pairs/"),
             data.table(form= paste0("diff~", mot[, paste0(V1, "+", V2)]), # Residuals / Motif pair w/0 interaction
                        save= paste0("predict_residuals__", mot[, paste0(V1, "+", V2)]),
                        folder= "db/linear_models/motif_pairs/"),
             data.table(form= paste0("log2FoldChange~median_L*median_R+", paste0(motMerge, collapse = "+")), # log2FC / Activity  and merged motif w/0 interaction
                        save= "predict_log2FoldChange__actL*R+mergedMot",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~median_L*median_R+", paste0(motMergeSom, collapse = "+")), # log2FC / Activity  and merged motifSOM w/0 interaction
                        save= "predict_log2FoldChange__actL*R+mergedMotSom",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("log2FoldChange~", paste0(motMerge, collapse = "+")), # log2FC / Merged motifs only w/0 interaction
                        save= "predict_log2FoldChange__mergedMot",
                        folder= "db/linear_models/global_models/"),
             #### Residuals prediction
             data.table(form= paste0("diff~", paste0(paste0(mot[, paste0(V1, "*", V2)], collapse= "+"))),  # Residuals / Motifs only w/ interaction
                        save= "predict_residuals__motL*R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("diff~", paste0(unique(unlist(mot)), collapse= "+")), # Residuals / Motifs only w/o interaction
                        save= "predict_residuals__motL+R",
                        folder= "db/linear_models/global_models/"),
             data.table(form= paste0("diff~", paste0(motMerge, collapse = "+")), # Residuals / Merged motifs only w/0 interaction
                        save= "predict_residuals__mergedMot",
                        folder= "db/linear_models/global_models/"))
mod <- mod[, .(group= unique(dat$active_plot_group)), (mod)]
mod[, save:= paste0(folder, gsub(": | ", "_", group), "__", save, ".rds"), .(folder, save, group)]

#-------------------------#
# MODELS
#-------------------------#
for(i in seq(nrow(mod)))
{
  if(!file.exists(mod[i, save]))
  {
    print(paste0(mod[i, save], " -->> START!"))
    .c <- lm(as.formula(mod[i, form]), dat[mod[i, group]])
    saveRDS(.c, mod[i, save])
    print(paste0(mod[i, save], " -->> DONE!"))
  }
}