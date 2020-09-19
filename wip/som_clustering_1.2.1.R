load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")
require(factoextra)
require(kohonen)
require(pheatmap)

som.model <- readRDS("Rdata/som_peSTARRSeq.rds")

# Add nodes to dat
cl <- data.table(uniq_ID= rownames(som.model$data$enh), node= som.model$unit.classif)
c_dat <- copy(dat)
c_dat[cl, node_L:= i.node, on= "enh_L==uniq_ID"]
c_dat[cl, node_R:= i.node, on= "enh_R==uniq_ID"]
c_dat <- c_dat[!is.na(node_L) & !is.na(node_R)]

# Residuals matrix 
mat <- dcast(c_dat, node_L~node_R, value.var = "diff", fun.aggregate = median, na.rm= T)

# Motifs
c_mot <- as.matrix(som.model$codes$motifs)
colnames(c_mot) <- sapply(colnames(c_mot), function(x) 
{
  cur <- strsplit(x, "_")[[1]]
  cur <- cur[length(cur)]
  return(ifelse(cur=="NA", gsub("_NA$", "", x), cur))
})
mot_ph <- pheatmap(c_mot, col= colorRampPalette(c("black", "blue", "yellow", "red"))(100), scale= "row", cutree_rows = 6)
mot_cl <- cutree(mot_ph$tree_row, 6)

# Features
c_feat <- as.matrix(som.model$codes$features)
feat_ph <- pheatmap(c_feat, col= colorRampPalette(c("black", "blue", "yellow", "red"))(100), scale= "row", cutree_rows = 4)
mot_cl <- cutree(mot_ph$tree_row, 6)


row_anot <- data.frame(individual_activity= som.model$codes$enh[, "median_L"], 
                       )

cl <- 



# Clustering
m1 <- as.matrix(dcast(c_dat, node_L~node_R, value.var = "diff", fun.aggregate= median, na.rm= T), 1)
m1[is.na(m1)] <- 0
clip <- quantile(m1, c(0.01, 0.99))
m1[m1<clip[1]] <- clip[1]
m1[m1>clip[2]] <- clip[2]

pdf("pdf_wip/heatmap_nodes.pdf", 8, 8)
par(mar= c(5,5,5,7))
r1 <- my_pheatmap(m1, cutree_rows = 7, cutree_cols = 7, lim= c(-2, 2))
dev.off()

# Clustering active pairs
sub <- c_dat[]
m2 <- as.matrix(dcast(sub, node_L~node_R, value.var = "diff", fun.aggregate= mean, na.rm= T), 1)
ranot <- as.matrix(som.model$codes$enh[, "median_L"])[as.numeric(rownames(m2)),, drop= F]


rownames(ranot) <- gsub("V", "", rownames(ranot))
colnames(ranot) <- "activity"
ranot <- as.data.frame(ranot)

pdf("pdf_wip/heatmap_nodes.pdf", 12, 12)

pheatmap(m2, annotation_row = ranot, annotation_col = ranot, 
         annotation_colors = list(activity= colorRampPalette(c("yellow", "blue"))(15)))


dev.off()







