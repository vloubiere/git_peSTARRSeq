setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(digest)
dir.create("db/linear_models/collapsed_data/", 
           showWarnings = F)
options(datatable.prettyprint.char=35L)

#--------------------------------#
# Import and compute extra variables
#--------------------------------#
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")
dat <- DT[!is.na(additive) & L!=R]
dat <- dat[(cdition== "vllib002" & group_L== "dev" & group_R== "dev" & median_L>1 & median_R>1) |
           (cdition== "vllib015" & group_L== "dev" & group_R== "dev") |
           (cdition== "vllib016" & group_L== "hk" & group_R== "hk")]

#--------------------------------#
# Build collapsed version
#--------------------------------#
cols <- c("cdition", "L", "median_L", grep("^motif.*_L$", names(dat), value = T))
.L <- dat[, .(side= "L", median= median(median_R), log2FoldChange= median(log2FoldChange), diff= median(diff)), cols]
setnames(.L, c("L", "median_L"), c("ID", "act"))
names(.L) <- gsub("(motif.*)_L", "\\1", names(.L))
cols <- c("cdition", "R", "median_R", grep("^motif.*_R$", names(dat), value = T))
.R <- dat[, .(side= "R", median= median(median_L), log2FoldChange= median(log2FoldChange), diff= median(diff)), cols]
setnames(.R, c("R", "median_R"), c("ID", "act"))
names(.R) <- gsub("(motif.*)_R", "\\1", names(.R))
coll <- rbind(.L, .R)

#--------------------------------#
# Motifs formulas
#--------------------------------#
# Using som code vectors
addSom <- paste0(grep("^motifSom__", names(coll), value = T), collapse= "+")

#-------------------------#
# Build models formulas
#-------------------------#
mod <- data.table(cdition= unique(coll$cdition))
mod <- mod[, .(dependent= c("log2FoldChange", "diff")), (mod)]
mod <- mod[, .(ind_variable1= c("act+median", "act*median", "act+median+side", "act*median*side")), (mod)]
mod <- mod[, .(ind_variable2= c(NA, addSom)), (mod)]
mod[, form:= paste0(dependent, "~", ind_variable1)]
mod[!is.na(ind_variable2), form:= paste0(form, "+", ind_variable2)]
mod[, file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/linear_models/collapsed_data/", 
                    cdition, "_", 
                    substr(form, 1, 50), "_", 
                    .GRP, ".rds"), (mod)]

#-------------------------#
# MODELS using full data
#-------------------------#
setkeyv(coll, "cdition")
for(i in seq(nrow(mod)))
{
  # Check if model already saved
  # file.remove(list.files("db/linear_models/collapsed_data/", full.names = T))
  if(!file.exists(mod[i,file]))
  {
    print(paste0(mod[i, file], " -->> START!"))
    # Build model
    .c <- lm(as.formula(mod[i, form]), coll[mod[i, cdition]])
    saveRDS(.c, mod[i,file])
    print(paste0(mod[i, file], " -->> DONE!"))
  }
}


mod[, c("R2", "equation", "PCC", "coeffs"):= {
  .c <- readRDS(file)
  .s <- summary(.c)
  .(.s$r.squared, 
    vl_model_equation(.c),
    cor(.c$model[,1], predict(.c)),
    .(as.data.table(.s$coefficients, keep.rownames= T)))
}, file]

saveRDS(mod, "Rdata/collapsed_data_modelling_final_table.rds")
