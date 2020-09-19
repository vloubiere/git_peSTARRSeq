load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")
source("/groups/stark/vloubiere/scripts/R_functions/my_plots.R")
require(pheatmap)
require(kohonen)

individual_act <- merge(unique(dat[!is.na(median_L), .(enh= enh_L, median_L)]), unique(dat[!is.na(median_R), .(enh= enh_R, median_R)]))
enh <- individual_act$enh
act_L <- dcast(dat[enh_L %in% enh & enh_R %in% enh], enh_L~enh_R, value.var= "log2FoldChange")
act_R <- dcast(dat[enh_L %in% enh & enh_R %in% enh], enh_R~enh_L, value.var= "log2FoldChange")
diff_L <- dcast(dat[enh_L %in% enh & enh_R %in% enh], enh_L~enh_R, value.var= "diff")
diff_R <- dcast(dat[enh_L %in% enh & enh_R %in% enh], enh_R~enh_L, value.var= "diff")
if(!exists("c_mot"))
{
  c_mot <- mot[!is.na(Dmel_prot) & gsub("A|C|G|T|R|Y|S|W|K|M|B|D|H|V|N", "", Dmel_prot) != "" & uniq_ID %in% individual_act$enh]
  c_mot[, check1:= length(which(low_motif_count>0))>25, motif]
  c_mot[, check2 := sum(low_motif_count, na.rm= T), motif]
  c_mot[, check3 := motif[which.max(check2)], Dmel_prot]
  c_mot <- c_mot[(check1) & motif==check3, !c("check1", "check2", "check3")]
  c_mot[, low_motif_count:= 
          {
            cur <- as.data.table(table(low_motif_count))
            cur[, check:= 1-cumsum(N)/sum(N)]
            if(any(cur$check<0.02))
            {
              sat <- as.numeric(cur[max(which(check>=0.02)), low_motif_count])
              low_motif_count[low_motif_count>=sat] <- sat
            }
            low_motif_count
          }, motif]
}
motif <- dcast(c_mot, uniq_ID~Dmel_prot, value.var= "low_motif_count")
cols <- setdiff(colnames(motif), "uniq_ID")
motif[, (cols) := lapply(.SD, function(x) log2(x+1)), .SDcols= cols]
cols <- c("uniq_ID", colnames(feat)[24:length(colnames(feat))])
features <- feat[uniq_ID %in% enh, ..cols]

obj <- list(individual_act= individual_act, act_L= act_L, act_R= act_R, diff_L= diff_L, diff_R= diff_R, motif= motif, features= features)
obj <- lapply(obj, function(x) 
{
  setkeyv(x, colnames(x)[1])
  return(as.matrix(x[individual_act$enh], 1))
})
sel_L <- apply(obj$diff_L, 1, function(x) length(which(!is.na(x)))>200)
sel_R <- apply(obj$diff_R, 1, function(x) length(which(!is.na(x)))>200)
sel <- which(sel_L & sel_R)
obj <- lapply(obj, function(x) x[sel,])

grid.size <- 12
mygrid <- somgrid(xdim = grid.size, ydim = grid.size, topo = 'hexagonal', toroidal = T)
set.seed(1)
if(!exists("som.model"))
{
  # som.model <- supersom(obj, whatmap = names(obj), user.weights = rep(1, length(obj)), grid = mygrid, maxNA.fraction = 0.99)
  som.model <- supersom(obj, whatmap = names(obj), user.weights = c(1,1,1,100,100,1,1), grid = mygrid, maxNA.fraction = 0.99)
}



Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))
par(mfrow= c(4,4))
plot(som.model, type = 'property', property = som.model$codes$individual_act[, "median_L"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$individual_act[, "median_R"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$act_L[,"dev_medium_A_00001"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$act_R[,"dev_medium_A_00001"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$diff_L[,"dev_medium_A_00001"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$diff_R[,"dev_medium_A_00001"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$motif[,"Trl"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$motif[,"cad"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$motif[,"Antennapedia"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[,"H3K27me3_WT_GSE41440"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[,"PH_WT_GSE60686"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[,"H3K27ac_dslacz_GSE81795"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[,"H3K4me1_dslacz_GSE81795"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[,"H3K4me2_dslacz_GSE81795"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[,"H3K4me3_dslacz_GSE81795"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[, "RNAPolII_WT_GSE41440"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[, "PHO_WT_GSE84502"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[, "GAF_untreated_GSE40646"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)
plot(som.model, type = 'property', property = som.model$codes$features[, "CTCF_control_GSE41354"], keepMargins = F, main = '', shape= "straight", palette.name= Cc)


c_dat <- copy(dat)















