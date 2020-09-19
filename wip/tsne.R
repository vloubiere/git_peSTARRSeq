# par(mfrow=c(3, 4))
# sapply(1:5, function(x) boxplot(diff~cl_R, c_dat[cl_L==x], notch= T))
# 
# c_dat[, effect_R:= median(diff, na.omit=T), .(cl_L, enh_R)]
# t1 <- unique(c_dat[, .(cl_L, enh_R, effect_R)])
# t2 <- as.matrix(dcast(t1[!is.na(cl_L)], cl_L~enh_R, value.var = "effect_R"), 1)
# 
# dev.off()
# par(mar= c(5,5,5,7))
# Cc <- c("black", "blue", "yellow")
# my_pheatmap(t(t2), cluster_cols= F, row_labels = F, clustering_method="ward.D", lim= c(0, 2.5), col= Cc, cutree_rows = 9)
# 
# 
# library(Rtsne)
# 
# test <- t(t2)
# set.seed(1)
# tsne <- Rtsne(test, dims = 2, perplexity= 25, verbose=TRUE, max_iter = 500)
# 
# ## Plotting
# par(mfrow= c(5, 5), las= 1)
# Cc <- colorRamp2(seq(quantile(test, 0.01), quantile(test, 0.99), length.out = 3), colors = c("black", "blue", "gold"))
# for(cl in colnames(test))
# {
#   plot(tsne$Y, main= paste0("left cluster ", cl), col= Cc(test[, cl]), pch= 19)
# }
# for(sel in sample(c_mot$motif, 18))
# {
#   cur <- log2(mot[motif==sel][match(rownames(test), uniq_ID), low_motif_count]+1)
#   Cc <- colorRamp2(c(0, max(cur)), colors = c("black", "gold"))
#   plot(tsne$Y, col= Cc(cur), pch= 19, main= mot[motif==sel, Dmel_prot][1])
# }
