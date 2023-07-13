setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import vllib002
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
cols <- names(dat)
# Add expected values and reorder columns
dat[, additive:= log2(2^indL+2^indR)]
dat[, multiplicative:= indL+indR]
setcolorder(dat,
            c(first(cols, n = length(cols)-2), "additive", "multiplicative", "predicted", "residuals"))
# Add actPairs expected values (smaller model)
actPairs <- readRDS("db/linear_models/FC_vllib002_actPairs_lm_predictions.rds")
dat[actPairs, predicted.actPairs:= i.predicted, on= c("L", "R")]
dat[actPairs, residuals.actPairs:= i.residuals, on= c("L", "R")]
# Add PCA on actPairs residuals
pca <- readRDS("db/pca/vllib002_actPairs_pca_residuals.rds")
dat[pca, c("PC1L", "PC1R"):= .(i.PC1L, i.PC1R), on= c("L", "R")]
# Add motif counts main motifs
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
sel <- c('AP-1'= "jaspar__MA0491.1",
         'Trl'= "flyfactorsurvey__Trl_FlyReg_FBgn0013263",
         'twist'= "flyfactorsurvey__twi_da_SANGER_5_FBgn0000413",
         'Dref'= "homer__AVYTATCGATAD_DREF")
mot <- vl_motif_counts(lib$oligo_full_sequence, 
                       sel = sel, 
                       genome = "dm3", 
                       collapse_overlapping = T)
setnames(mot, names(sel))
motL <- mot[match(dat$L, lib$ID_vl),]
setnames(motL, function(x) paste0(x, "_L"))
dat <- cbind(dat, motL)
motR <- mot[match(dat$R, lib$ID_vl),]
setnames(motR, function(x) paste0(x, "_R"))
dat <- cbind(dat, motR)
saveRDS(dat,
        "Rdata/vlli002_final_table.rds")
