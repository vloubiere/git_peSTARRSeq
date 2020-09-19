load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")
require(factoextra)
require(kohonen)

if(!exists("c_mot"))
{
  # c_mot <- mot[!is.na(Dmel_prot) & gsub("A|C|G|T|R|Y|S|W|K|M|B|D|H|V|N", "", Dmel_prot) != ""]
  c_mot <- mot[(!is.na(ID_vl) & group=="dev") | BA_group=="NegativeRegions"]
  c_mot[is.na(group), group:= "control"]
  # Dev TFs
  c_mot[, c("ORdev", "pvaldev"):= fisher.test(table(group=="dev", low_motif_count>0))[c("estimate", "p.value")], motif]
  padj <- unique(c_mot[, .(motif, pvaldev)])[, .(motif, padj= p.adjust(pvaldev))]
  c_mot[padj, padjdev:= i.padj, on= "motif"]
  
  # Selection
  c_mot <- c_mot[(ORdev>1 & padjdev<0.05)]
  
  # Final filters
  c_mot[, check1 := sum(low_motif_count, na.rm= T), motif]
  c_mot[, check2 := motif[which.max(check1)], .(Dmel_prot, motif)]
  c_mot <- c_mot[motif==check2, !c("check1", "check2")]
  c_mot <- dcast(c_mot, uniq_ID~motif+Dmel_prot, value.var = "low_motif_count", fun.aggregate = function(x) log2(x+1), fill= 0)
  setkey(c_mot, uniq_ID)
}

if(!exists("c_feat"))
{
  cols <- c("uniq_ID", grep("GSE", colnames(feat), value= T))
  c_feat <- feat[!is.na(ID_vl)]
  c_feat <- c_feat[, ..cols]
  colnames(c_feat)[-1] <- sapply(colnames(c_feat)[-1], function(x) strsplit(x, "_")[[1]][1])
  setkey(c_feat, uniq_ID)
}

c_dat <- dat[act_group=="active~active"]

mat <- as.matrix(dcast(c_dat, enh_L~enh_R, value.var= "diff"), 1)
mat <- mat[apply(mat, 1, function(x) length(na.omit(x))>0.6*length(x)),]
mat <- mat[, apply(mat, 2, function(x) length(na.omit(x))>0.6*length(x))]
lim <- quantile(mat, c(0.01, 0.99), na.rm= T)
mat[mat < lim[1]] <- lim[1]
mat[mat > lim[2]] <- lim[2]

pdf("pdf_wip/heatmap_peSTARR_residuals_with_motifs.pdf", 24, 21)
layout(matrix(c(3,4,1,2), byrow = T, ncol= 2), widths = c(1, 0.7), heights = c(0.7, 1))
par(mar= c(7, 7, 2, 6))

res <- my_pheatmap(mat, clustering_dist_rows = "maximum", clustering_dist_cols = "maximum",
                   clustering_method_cols = "ward.D2", clustering_method_rows = "ward.D2", 
                   cluster_cols = 5, cutree_rows = 7, cutree_cols = 7, lwd= 2)

c_mot_L <- as.matrix(c_mot[unique(res[, .(row_labels, row_order)])[order(row_order), row_labels]], 1)
par(mar= c(7, 1, 2, 3))
my_pheatmap(c_mot_L, cluster_rows= F, row_labels = F, scale= "column", cutree_cols = 6, col_labels = T, lwd= 1)

c_mot_R <- t(as.matrix(c_mot[unique(res[, .(col_labels, col_order)])[order(col_order), col_labels]], 1))
par(mar= c(1, 7, 2, 6))
my_pheatmap(c_mot_R, cluster_cols= F, row_labels = T, scale= "row", cutree_rows = 6, col_labels = F, lwd= 1)
dev.off()

pdf("pdf_wip/heatmap_peSTARR_residuals_with_features.pdf", 24, 21)
layout(matrix(c(3,4,1,2), byrow = T, ncol= 2), widths = c(1, 0.7), heights = c(0.7, 1))
par(mar= c(7, 7, 2, 6))

res <- my_pheatmap(mat, clustering_dist_rows = "maximum", clustering_dist_cols = "maximum",
                   clustering_method_cols = "ward.D2", clustering_method_rows = "ward.D2", 
                   cluster_cols = 5, cutree_rows = 7, cutree_cols = 7, lwd= 2)

c_feat_L <- as.matrix(c_feat[unique(res[, .(row_labels, row_order)])[order(row_order), row_labels]], 1)
par(mar= c(7, 1, 2, 3))
my_pheatmap(c_feat_L, cluster_rows= F, row_labels = F, cutree_cols = 6, col_labels = T, lwd= 1, scale= "column")

c_feat_R <- t(as.matrix(c_feat[unique(res[, .(col_labels, col_order)])[order(col_order), col_labels]], 1))
par(mar= c(1, 7, 2, 6))
my_pheatmap(c_feat_R, cluster_cols= F, row_labels = T, cutree_rows = 6, col_labels = F, lwd= 1, scale= "row")
dev.off()



pdf("pdf_wip/heatmap_peSTARR_residuals_with_motifs_2.pdf", 24, 21)
layout(matrix(c(2,4,3,1), byrow = T, ncol= 2), widths = c(1, 0.7), heights = c(0.7, 1))

c_mot_L <- as.matrix(c_mot[unique(res[, .(row_labels, row_order)])[order(row_order), row_labels]], 1)
par(mar= c(7, 1, 2, 3))
res_L <- my_pheatmap(c_mot_L, row_labels = F, scale= "column", cutree_cols = 8, col_labels = T, lwd= 1)

c_mot_R <- t(as.matrix(c_mot[unique(res[, .(col_labels, col_order)])[order(col_order), col_labels]], 1))
par(mar= c(1, 7, 2, 6))
res_R <- my_pheatmap(c_mot_R, row_labels = T, scale= "row", cutree_rows = 8, col_labels = F, lwd= 1)

c_mat <- as.matrix(dcast(dat, enh_L~enh_R, value.var = "diff"), 1)
# c_mat <- as.matrix(dcast(dat, enh_L~enh_R, value.var = "log2FoldChange"), 1)
c_mat <- c_mat[match(res_L[, row_labels[1], keyby= row_order]$V1, rownames(c_mat)),]
c_mat <- c_mat[, match(res_R[, col_labels[1], keyby= col_order]$V1, colnames(c_mat))]
par(mar= c(7, 7, 2, 6))
my_pheatmap(c_mat, cluster_rows = F, cluster_cols = F)
dev.off()

pdf("pdf_wip/heatmap_peSTARR_residuals_with_features_2.pdf", 24, 21)
layout(matrix(c(2,4,3,1), byrow = T, ncol= 2), widths = c(1, 0.7), heights = c(0.7, 1))

c_feat_L <- as.matrix(c_feat[unique(res[, .(row_labels, row_order)])[order(row_order), row_labels]], 1)
c_feat_L[is.na(c_feat_L)] <- 0
par(mar= c(7, 1, 2, 3))
res_L <- my_pheatmap(c_feat_L, row_labels = F, scale= "column", cutree_cols = 8, col_labels = T, lwd= 1)

c_feat_R <- t(as.matrix(c_feat[unique(res[, .(col_labels, col_order)])[order(col_order), col_labels]], 1))
c_feat_R[is.na(c_feat_R)] <- 0
par(mar= c(1, 7, 2, 6))
res_R <- my_pheatmap(c_feat_R, row_labels = T, scale= "row", cutree_rows = 8, col_labels = F, lwd= 1)

# c_mat <- as.matrix(dcast(dat, enh_L~enh_R, value.var = "diff"), 1)
c_mat <- as.matrix(dcast(dat, enh_L~enh_R, value.var = "log2FoldChange"), 1)
c_mat <- c_mat[match(res_L[, row_labels[1], keyby= row_order]$V1, rownames(c_mat)),]
c_mat <- c_mat[, match(res_R[, col_labels[1], keyby= col_order]$V1, colnames(c_mat))]
par(mar= c(7, 7, 2, 6))
my_pheatmap(c_mat, cluster_rows = F, cluster_cols = F)
dev.off()









