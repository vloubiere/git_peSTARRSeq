load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")
require(factoextra)
require(kohonen)

if(!exists("c_mot"))
{
  # c_mot <- mot[!is.na(Dmel_prot) & gsub("A|C|G|T|R|Y|S|W|K|M|B|D|H|V|N", "", Dmel_prot) != ""]
  c_mot <- mot[!is.na(ID_vl) | BA_group=="NegativeRegions"]
  c_mot[is.na(group), group:= "control"]
  # Dev TFs
  c_mot[, c("ORdev", "pvaldev"):= fisher.test(table(group=="dev", low_motif_count>0))[c("estimate", "p.value")], motif]
  padj <- unique(c_mot[, .(motif, pvaldev)])[, .(motif, padj= p.adjust(pvaldev))]
  c_mot[padj, padjdev:= i.padj, on= "motif"]
  # hk
  c_mot[, c("ORhk", "pvalhk"):= fisher.test(table(group=="hk", low_motif_count>0))[c("estimate", "p.value")], motif]
  padj <- unique(c_mot[, .(motif, pvalhk)])[, .(motif, padj= p.adjust(pvalhk))]
  c_mot[padj, padjhk:= i.padj, on= "motif"]
  # heatshock
  c_mot[, c("OROSC", "pvalOSC"):= fisher.test(table(group=="OSC", low_motif_count>0))[c("estimate", "p.value")], motif]
  padj <- unique(c_mot[, .(motif, pvalOSC)])[, .(motif, padj= p.adjust(pvalOSC))]
  c_mot[padj, padjOSC:= i.padj, on= "motif"]
  # ecdysone
  c_mot[, c("ORecdysone", "pvalecdysone"):= fisher.test(table(group=="ecdysone", low_motif_count>0))[c("estimate", "p.value")], motif]
  padj <- unique(c_mot[, .(motif, pvalecdysone)])[, .(motif, padj= p.adjust(pvalecdysone))]
  c_mot[padj, padjecdysone:= i.padj, on= "motif"]
  # OSC
  c_mot[, c("ORheatshock", "pvalheatshock"):= fisher.test(table(group=="heatshock", low_motif_count>0))[c("estimate", "p.value")], motif]
  padj <- unique(c_mot[, .(motif, pvalheatshock)])[, .(motif, padj= p.adjust(pvalheatshock))]
  c_mot[padj, padjheatshock:= i.padj, on= "motif"]
  
  # Selection
  c_mot <- c_mot[(ORdev>1 & padjdev<0.001) |
                 (ORhk>1 & padjhk<0.001) |
                 (OROSC>1 & padjOSC<0.001) |
                 (ORecdysone>1 & padjecdysone<0.001) |
                 (ORheatshock>1 & padjheatshock<0.001)]
  
  # Final filters
  c_mot[, check1 := sum(low_motif_count, na.rm= T), motif]
  c_mot[, check2 := motif[which.max(check1)], Dmel_prot]
  c_mot <- c_mot[motif==check2, !c("check1", "check2")]
  c_mot <- dcast(c_mot, uniq_ID~motif+Dmel_prot, value.var = "low_motif_count", fun.aggregate = function(x) log2(x+1), fill= 0)
  colnames(c_mot) <- gsub(" ", ".", colnames(c_mot))
  setkey(c_mot, uniq_ID)
}

if(!exists("c_feat"))
{
  c_feat <- feat[!is.na(ID_vl)]
  cols <- c("uniq_ID", grep("GSE", colnames(feat), value= T))
  cols <- grep("H3K4me3_WT_GSE41440|Trr_WT_GSE41440|PH_ML", cols, invert = T, value = T) # problem with this sample / irrelevant
  c_feat <- c_feat[, ..cols]
  colnames(c_feat)[-1] <- sapply(colnames(c_feat)[-1], function(x) strsplit(x, "_")[[1]][1])
  setkey(c_feat, uniq_ID)
}

if(!exists("obj"))
{
  clip <- function(mat, qmin= 0.001, qmax= 0.999)
  {
    .min <- quantile(mat, qmin, na.rm= T)
    .max <- max(mat, qmax, na.rm= T)
    mat[mat < .min] <- .min
    mat[mat > .max] <- .max
    return(mat)
  }
  enh <- data.table(uniq_ID= sort(unique(unlist(c(dat[!is.na(median_L) & !is.na(median_R), .(enh_L, enh_R)])))))
  c_dat <- dat[enh_L %in% enh$uniq_ID & enh_R %in% enh$uniq_ID & !is.na(median_L) & ! is.na(median_R)]
  enh[, median_L:= c_dat[.BY, unique(median_L), on= "enh_L==uniq_ID"], uniq_ID]
  enh[, median_R:= c_dat[.BY, unique(median_R), on= "enh_R==uniq_ID"], uniq_ID]
  
  class <- as.factor(sapply(enh$uniq_ID, function(x) strsplit(x, "_")[[1]][1]))
  motifs <- as.matrix(c_mot[enh[, uniq_ID], , on = "uniq_ID"], 1)
  features <- as.matrix(c_feat[enh[, uniq_ID], , on = "uniq_ID"], 1)
  
  act_L <- as.matrix(dcast(c_dat, enh_L~enh_R, value.var = "log2FoldChange"), 1)
  act_L <- act_L[match(enh$uniq_ID, rownames(act_L)),]
  act_L <- act_L[,match(enh$uniq_ID, colnames(act_L))]
  act_L <- clip(act_L)
  
  act_R <- as.matrix(dcast(c_dat, enh_R~enh_L, value.var = "log2FoldChange"), 1)
  act_R <- act_R[match(enh$uniq_ID, rownames(act_R)),]
  act_R <- act_R[,match(enh$uniq_ID, colnames(act_R))]
  act_R <- clip(act_R)
  
  diff_L <- as.matrix(dcast(c_dat, enh_L~enh_R, value.var = "diff"), 1)
  diff_L <- diff_L[match(enh$uniq_ID, rownames(diff_L)),]
  diff_L <- diff_L[,match(enh$uniq_ID, colnames(diff_L))]
  diff_L <- clip(diff_L)
  
  diff_R <- as.matrix(dcast(c_dat, enh_R~enh_L, value.var = "diff"), 1)
  diff_R <- diff_R[match(enh$uniq_ID, rownames(diff_R)),]
  diff_R <- diff_R[,match(enh$uniq_ID, colnames(diff_R))]
  diff_R <- clip(diff_R)
  
  cor_L <- cor(as.matrix(dcast(c_dat, enh_R~enh_L, value.var = "diff"), 1), use= "pairwise.complete.obs")
  cor_L <- cor_L[match(enh$uniq_ID, rownames(cor_L)),]
  cor_L <- cor_L[,match(enh$uniq_ID, colnames(cor_L))]
  
  cor_R <- cor(as.matrix(dcast(c_dat, enh_L~enh_R, value.var = "diff"), 1), use= "pairwise.complete.obs")
  cor_R <- cor_R[match(enh$uniq_ID, rownames(cor_R)),]
  cor_R <- cor_R[,match(enh$uniq_ID, colnames(cor_R))]
  
  obj <- list(enh= as.matrix(enh, 1), class= class, motifs= motifs, features= features,
              act_L= act_L, act_R= act_R, diff_L= diff_L, diff_R= diff_R,
              cor_L= cor_L, cor_R= cor_R)
}

if(!file.exists("Rdata/som_peSTARRSeq.rds"))
{
  mygrid <- somgrid(xdim= 6, ydim= 6, topo = 'hexagonal', toroidal = T)
  set.seed(123)
  som.model <- supersom(obj, whatmap = names(obj), user.weights = c(1,1,1,1,1,1,1,1,8,8), grid = mygrid, maxNA.fraction = 0.9, rlen = 1000)
  
  ## use hierarchical clustering to cluster the codebook vectors
  fviz_nbclust(as.matrix(object.distances(som.model, "codes")), FUN = hcut, method = "wss", k.max = 30)
  som.model$hclust <- cutree(hclust(object.distances(som.model, "codes")), 9)
  
  # save
  saveRDS(som.model, "Rdata/som_peSTARRSeq.rds")
}

# diagnostic plots
som.model <- readRDS("Rdata/som_peSTARRSeq.rds")
Cc <- colorRampPalette(c("black", "blue", "yellow"))

pdf("pdf_wip/som_diagnostic.pdf", 15, 15)
par(mfrow= c(4, 4))
for(type in c("changes", "dist.neighbours", "counts", "mapping", "quality"))
{
  plot(som.model, type= type, palette.name = Cc, shape= "straight", border= "black")
  if(type != "changes")
  {
    add.cluster.boundaries(som.model, som.model$hclust, lwd= 3, col= "red")
  }
}
for(what in names(som.model$codes))
{
  plot(som.model, type= "codes", whatmap= what, shape= "straight", border= "black")
  add.cluster.boundaries(som.model, som.model$hclust, lwd= 3, col= "red")
}
dev.off()
