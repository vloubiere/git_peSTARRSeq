setwd("/groups/stark/vloubiere/projects/pe_STARRSeq")
require(data.table)
require(motifmatchr)
require(TFBSTools)
require(kohonen)
dir_pdf <- "pdf/motifs_som_clustering"
dir.create(dir_pdf, showW)

#--------------------------------------#
# Motif counts
#--------------------------------------#
lib <- readRDS("Rdata/uniq_300bp_enhancers.rds")
lib <- lib[!detail=="ecoli"]

# Count hits low cutoff
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
sel <- as.data.table(TF_clusters_PWMs$metadata)[X..motif_collection_name %in% c("jaspar", "bergman", "cisbp", "hocomoco", "homer", "flyfactorsurvey"), motif_name]
sel <- name(TF_clusters_PWMs$All_pwms_log_odds) %in% sel
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], 
                   GRanges(lib$seqnames, IRanges(lib$start, lib$end)), 
                   genome= "dm3", 
                   p.cutoff= 5e-4, 
                   bg="even", 
                   out= "scores")
counts <- as.matrix(motifCounts(hit))
colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds)[sel]
rownames(counts) <- lib$ID
counts <- as.data.table(counts, keep.rownames = T)

#--------------------------------------#
# Select motifs enriched in at least one of the lib groups
#--------------------------------------#
dat <- melt(counts, 
            id.vars= "rn")
dat[lib, group:= i.group, on= "rn==ID"]

cat_counts <- dat[, .(no_motif= length(which(value==0)), all= .N), .(variable, group)]
cat_counts[, motif:= all-no_motif]
fisher <- merge(cat_counts[group!="control"], 
                cat_counts[group=="control", !"group"], 
                by= "variable",
                suffixes= c("_gp", "_control"))
fisher[, c("estimate", "pvalue"):= fisher.test(matrix(c(no_motif_control,
                                                        motif_control,
                                                        no_motif_gp, 
                                                        motif_gp), 
                                                      ncol= 2, 
                                                      byrow = T))[c("estimate", "p.value")], .(variable, group)]
fisher[, padj:= p.adjust(pvalue, method= "fdr")]
sel_motifs <- unique(fisher[padj<0.00001 & estimate>2, variable])

#--------------------------------------#
# SOM clustering
#--------------------------------------#
dmat <- dcast(dat[variable %in% sel_motifs], 
              variable~rn, 
              value.var = "value")
mat <- log2(as.matrix(dmat, 1)+1)
mat <- scale(mat) 
mat <- apply(mat, 2, function(x)
{
  lim <- quantile(x, 0.01, 0.99)
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  return(x)
})
mygrid <- somgrid(xdim= 4, 
                  ydim= 4, 
                  topo = 'hexagonal', 
                  toroidal = T)
set.seed(1234)
som.model <- supersom(mat, 
                      grid = mygrid, 
                      rlen = 500)

# Diag plots
pdf(paste0(dir_pdf, "/SOM_motifs_clustering_diag_plots.pdf"))
par(mfrow= c(4, 4))
for(type in c("changes", "dist.neighbours", "counts", "mapping", "quality"))
  plot(som.model, 
       type= type, 
       palette.name = colorRampPalette(c("blue", "yellow")), 
       shape= "straight", 
       border= "black")
dev.off()

# Heatmap
pdf(paste0(dir_pdf, "/SOM_motifs_clustering_heatmap.pdf"))
pl <- mat[order(som.model$unit.classif),]
vl_heatmap(pl, 
           cluster_rows = F, 
           col = colorRampPalette(c("blue", "yellow", "white"))(100))
abline(h= 1-which(diff(som.model$unit.classif[order(som.model$unit.classif)])!=0)/nrow(mat), lwd= 1, col= "white")
dev.off()

# Identify representative motifs
cl <- data.table(cl= som.model$unit.classif, 
                 dist= som.model$distances, 
                 motif= rownames(som.model$data[[1]]))
cl[as.data.table(TF_clusters_PWMs$metadata), Dmel:= i.Dmel, on= "motif==motif_name"]
cl[as.data.table(TF_clusters_PWMs$metadata), BA_cluster:= i.Motif_cluster_name, on= "motif==motif_name"]

# SAVE
som.model$info <- cl
saveRDS(som.model, "Rdata/som_clustering_motifs_300bp_enhancers.rds")

