setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")

# Dcast and clean
mat <- dcast(dat, -indL+L~indR+R, value.var = "residuals", sep = "__")
mat <- as.matrix(mat[, -1], 1)
colnames(mat) <- unlist(tstrsplit(colnames(mat), "__", keep= 2))
while(sum(is.na(mat))>0.05*nrow(mat)*ncol(mat))
{
  mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
  mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
}
# Complete NAs before left PCA
matL <- missMDA::imputePCA(mat, 
                           ncp = missMDA::estim_ncpPCA(mat)$ncp, 
                           scale = TRUE)
pcaL <- as.data.table(prcomp(matL[[1]])$x[, "PC1", drop= F], keep.rownames = "name")
# Complete NAs before Right PCA
matR <- missMDA::imputePCA(t(mat), 
                           ncp = missMDA::estim_ncpPCA(t(mat))$ncp, 
                           scale = TRUE)
pcaR <- as.data.table(prcomp(matR[[1]])$x[, "PC1", drop= F], keep.rownames = "name")

# Save object
pca <- CJ(L= pcaL$name,
          R= pcaR$name,
          unique = T)
pca[pcaL, PC1L:= i.PC1, on= "L==name"]
pca[pcaR, PC1R:= i.PC1, on= "R==name"]
saveRDS(pca, "db/pca/vllib002_pca_residuals.rds")
