setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(digest)
dir.create("db/linear_models/full_data/", 
           showWarnings = F)
options(datatable.prettyprint.char=50L)

#--------------------------------#
# Import and compute extra variables
#--------------------------------#
if(!exists("DT"))
  DT <- readRDS("Rdata/final_results_table.rds")
dat <- DT[!is.na(additive) & L!=R]
dat <- dat[(cdition== "vllib002" & group_L== "dev" & group_R== "dev" & median_L>1 & median_R>1) |
           (cdition== "vllib015" & group_L== "dev" & group_R== "dev" & median_L>1 & median_R>1) |
           (cdition== "vllib016" & group_L== "hk" & group_R== "hk" & median_L>1 & median_R>1)]
# Sum left right motifs
mot <- data.table(motif_LR= paste0(unique(gsub("_L|_R", "", grep("motif", names(dat), value = T))), "_LR"))
mot[, motL:= gsub("R$", "", motif_LR)]
mot[, motR:= gsub("LR$", "R", motif_LR)]
mot[, {
  dat[, (motif_LR):= get(motL)+get(motR)]
  print("")
}, (mot)]
# Distance between enhancers
dat[, dist:= abs(start(resize(GRanges(dat$coor_L), 1, "center"))-start(resize(GRanges(dat$coor_R), 1, "center")))]
dat[unlist(tstrsplit(dat$coor_L, ':', keep= 1))!=unlist(tstrsplit(dat$coor_R, ':', keep= 1)), dist:= max(dat$dist)]
dat[, dist:= log10(dist+1)]

#--------------------------------#
# Motifs formulas
#--------------------------------#
# Using som code vectors
som <- data.table(motL= grep("^motifSom__.*_L", names(dat), value = T),
                   motR= grep("^motifSom__.*_R", names(dat), value = T))
addSom <- paste0(unique(unlist(som)), collapse = "+")
multSom <- paste0(som[, apply(.SD, 1, paste0, collapse= "*")], collapse= "+")
# Using top motifs
top <- data.table(motL= grep("^motif__.*_L", names(dat), value = T),
                  motR= grep("^motif__.*_R", names(dat), value = T))
addTop <- paste0(unique(unlist(top)), collapse = "+")
multTop <- paste0(top[, apply(.SD, 1, paste0, collapse= "*")], collapse= "+")
# Using summed motifs
addSomLR <- paste0(grep("^motifSom__.*_LR", names(dat), value = T), collapse = "+")
addTopLR <- paste0(grep("^motif__.*_LR", names(dat), value = T), collapse = "+")

#-------------------------#
# Build models formulas
#-------------------------#
mod <- data.table(cdition= unique(dat$cdition))
mod <- mod[, .(dependent= c("log2FoldChange", "diff")), (mod)]
mod <- mod[, .(ind_variable1= c("median_L+median_R",
                                "median_L+median_R+dist",
                                "median_L*median_R",
                                "median_L*median_R+dist")), (mod)]
mod <- mod[, .(ind_variable2= c(NA,
                                addSom,
                                multSom,
                                addTop,
                                multTop,
                                addSomLR,
                                addTopLR)), (mod)]
mod[, form:= paste0(dependent, "~", ind_variable1)]
mod[!is.na(ind_variable2), form:= paste0(form, "+", ind_variable2)]
mod[, file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/linear_models/full_data/",
                    cdition, "_", 
                    substr(form, 1, 50), "_", 
                    .GRP, ".rds"), (mod)]


#-------------------------#
# MODELS using full data
#-------------------------#
setkeyv(dat, "cdition")
for(i in seq(nrow(mod)))
{
  # Check if model already saved
  # file.remove(list.files("db/linear_models/full_data/", full.names = T))
  if(!file.exists(mod[i,file]))
  {
    print(paste0(mod[i, file], " -->> START!"))
    # Build model
    .c <- lm(as.formula(mod[i, form]), dat[mod[i, cdition]])
    #Save
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

saveRDS(mod, "Rdata/full_data_modelling_final_table.rds")
