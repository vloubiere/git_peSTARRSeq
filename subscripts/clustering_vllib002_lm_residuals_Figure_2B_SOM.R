setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(kohonen)

# Import data
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh."]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")$model
dat[, diff:= log2FoldChange-predict(model, new= dat)]

# SOM clustering
trainL <- dcast(dat,
                L~R, 
                value.var = "diff")
trainL <- as.matrix(trainL, 1)
somGrid <- somgrid(xdim = 10, ydim = 10, topo= 'hexagonal', toroidal = T)
set.seed(1)
init <- sample(nrow(trainL), somGrid$xdim*somGrid$ydim)
init <- trainL[init, , drop=F]
somL <- supersom(data= trainL, 
                 grid= somGrid, 
                 init= init, 
                 rlen= 100, 
                 keep.data = TRUE,
                 maxNA.fraction = .2)

trainR <- dcast(dat,
                R~L, 
                value.var = "diff")
trainR <- as.matrix(trainR, 1)
somGrid <- somgrid(xdim = 10, ydim = 10, topo= 'hexagonal', toroidal = T)
set.seed(1)
init <- sample(nrow(trainR), somGrid$xdim*somGrid$ydim)
init <- trainR[init, , drop=F]
somR <- supersom(data= trainR, 
                 grid= somGrid, 
                 init= init, 
                 rlen= 100, 
                 keep.data = TRUE,
                 maxNA.fraction = .2)

# Codes clustering
set.seed(1)
clL <- data.table(name= rownames(somL$data[[1]]),
                  cl= kmeans(somL$codes[[1]], 4)$cluster[somL$unit.classif])
clL[, ord:= median(somL$data[[1]][rownames(somL$data[[1]]) %in% name,], na.rm = T), cl]
clL[, cl:= .GRP, keyby= ord]
clR <- data.table(name= rownames(somR$data[[1]]),
                  cl= kmeans(somR$codes[[1]], 4)$cluster[somR$unit.classif])
clR[, ord:= median(somR$data[[1]][rownames(somR$data[[1]]) %in% name,], na.rm = T), cl]
clR[, cl:= .GRP, keyby= ord]

# Heatmap
mat <- trainL[match(clL$name, rownames(trainL)), match(clR$name, colnames(trainL))]
cl <- vl_heatmap(mat,
                 row_clusters = clL$cl,
                 col_clusters = clR$cl,
                 cluster_rows = F, 
                 cluster_cols = F,
                 breaks = c(-2,-0.25,0.25,2), 
                 col= c("cornflowerblue", "white", "white", "tomato"), 
                 legend_title = "Obs./Exp. (log2)", 
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
plot(cl)
res <- merge(cl$rows, cl$cols, by= "name")
tab <- matrix(c(table(res$cl.x, res$cl.y)), ncol= 4)
.t <- chisq.test(tab)
set.seed(1)
rdm <- matrix(c(table(sample(res$cl.x), sample(res$cl.y))), ncol= 4)
OR <- log2((tab+1)/(rdm+1))
pl <- vl_heatmap(x = .t$residuals,
                 breaks = c(-3,0,3),
                 cluster_rows = F, 
                 cluster_cols = F, 
                 legend_title = "OR (log2)")
text(rep(1:4, each= 4), rep(4:1, 4), unlist(tab))


# saveRDS(cl, 
#         "Rdata/vllib002_clustering_expected_scores_draft_figure_pca.rds")