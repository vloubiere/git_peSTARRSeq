load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")
source("/groups/stark/vloubiere/scripts/R_functions/my_plots.R")
require(pheatmap)
require(gridExtra)

#### Select subset of active enhancers
c_dat <- unique(dat[act_group=="active~active", enh_L])
c_dat <- c_dat[c_dat %in% unique(dat[act_group=="active~active", enh_R])]
c_dat <- dat[enh_L %in% c_dat & enh_R %in% c_dat]

##### Clustering based on motifs
c_mot <- mot[!is.na(Dmel_prot) & gsub("A|C|G|T|R|Y|S|W|K|M|B|D|H|V|N", "", Dmel_prot) != "" & uniq_ID %in% unique(c(c_dat$enh_L,c_dat$enh_R))]
c_mot[, check1:= length(which(low_motif_count>0))>50, motif]
c_mot[, check2 := sum(low_motif_count, na.rm= T), motif]
c_mot[, check3 := motif[which.max(check2)], Dmel_prot]
c_mot <- c_mot[(check1) & motif==check3, !c("check1", "check2", "check3")]

mat <- as.matrix(dcast(c_mot, uniq_ID~Dmel_prot, value.var = "low_motif_count"), 1)
set.seed(12345)
km <- pheatmap(mat, kmeans_k = 6, scale= "column", col= colorRampPalette(c("cornflowerblue", "white", "tomato"))(100), 
               breaks= seq(-3, 3, length.out = 101), filename = "pdf_main/motif_clustering_peSTARRSeq_1.0.pdf", height = 6, width = 30)

cl <- data.table(uniq_ID= names(km$kmeans$cluster), cl= km$kmeans$cluster)

# Add to c_mot & c_dat
c_mot[cl, uniq_cl:= i.cl, on= "uniq_ID"]
c_dat[cl, cl_L:= i.cl, on= "enh_L==uniq_ID"]
c_dat[cl, cl_R:= i.cl, on= "enh_R==uniq_ID"]
c_dat[, uniq_cl:= paste0(cl_L, "~", cl_R)]

#### cl~cl additivity
pl <- copy(c_dat)
pl[, ord:= median(diff), uniq_cl]
pl <- pl[order(ord)]
pl[, ord:= .GRP+1, uniq_cl]
all <- copy(c_dat)
all[, c("uniq_cl", "ord"):= .("all", 1)]
pl <- rbind(all, pl)
pl[, c("cl_L", "cl_R"):= lapply(.(cl_L, cl_R), as.character)]
pl[uniq_cl=="all", cl_L:= "all"]
pl[uniq_cl=="all", cl_R:= "all"]

pdf("pdf_main/motifs_clusters_additivity_boxplots.pdf", width= 10, height= 12)
layout(matrix(c(1,2,3,3), ncol= 2, byrow= T))
par(las= 1)
boxplot(diff~cl_L, pl, ylim= c(-0.25, 6), notch= T, xlab= "cluster L", ylab= "observed vs expected additive", outline= F)
my_pval_plot(diff~cl_L, pl, pairs = lapply(1:6, function(x) c(x, 7)), y_adj = seq(1, 20, length.out = 6))
abline(h= median(pl[cl_L=="all", diff]), lty= 2)

boxplot(diff~cl_R, pl, ylim= c(-0.25, 6), notch= T, xlab= "cluster R", ylab= "observed vs expected additive", outline= F)
my_pval_plot(diff~cl_R, pl, pairs = lapply(1:6, function(x) c(x, 7)), y_adj = seq(1, 20, length.out = 6))
abline(h= median(pl[cl_R=="all", diff]), lty= 2)

box <- boxplot(diff~ord, pl, notch= T, outline= F, names= unique(pl$uniq_cl), las= 2, xlab= "L~R clusters cobinations", ylab= "observed vs expected additive")
abline(h= median(pl[uniq_cl=="all", diff]), lty= 2)
pval <- pl[uniq_cl!= "all", wilcox.test(diff, pl[uniq_cl=="all", diff])$p.value, uniq_cl]$V1
text(seq(2, ncol(box$stats)), box$stats[5,][-1], my_pval_format(pval)$p.value, cex= 0.5*my_pval_format(pval)$cex, pos= 3, offset= 0.2)
dev.off()


#
if(!exists("mot_cmbn"))
{
  # mot_cmbn <- CJ(mot_L= colnames(km$kmeans$centers)[km$kmeans$centers["1",]>0.5], mot_R= colnames(km$kmeans$centers)[km$kmeans$centers["1",]>0.5], unique = T)
  mot_cmbn <- CJ(mot_L= colnames(km$kmeans$centers)[apply(km$kmeans$centers,2, function(x) any(x>1))], mot_R= colnames(km$kmeans$centers)[apply(km$kmeans$centers,2, function(x) any(x>0.5))], unique = T)
  enh_cmbn <- unique(c_dat[, .(enh_L, enh_R)])
  mot_cmbn <- mot_cmbn[, enh_cmbn[, .(enh_L, enh_R)], .(mot_L, mot_R)]
  mot_cmbn[c_mot, count_L:= i.low_motif_count, on= c("enh_L==uniq_ID", "mot_L==Dmel_prot")]
  mot_cmbn[c_mot, count_R:= i.low_motif_count, on= c("enh_R==uniq_ID", "mot_R==Dmel_prot")]
  mot_cmbn[c_dat, diff:= i.diff, on= c("enh_L", "enh_R")]
  mot_cmbn[, c("tval", "pval") := 
             {
               res <- lm(diff~count_L*count_R)
               list(summary(res)$coefficients[4, 3], summary(res)$coefficients[4, 4])
             }, .(mot_L, mot_R)]
}
res <- unique(mot_cmbn[, .(mot_L, mot_R, tval, pval)])




















