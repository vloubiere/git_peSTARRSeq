setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(kohonen)

# Import data
# dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh."]
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
pred <- readRDS("Rdata/CV_linear_model_vllib002.rds")$pred
dat[pred, residuals:= log2FoldChange-i.predicted, on= c("L", "R")]

# SOM clustering
trainL <- dcast(dat,
                L~R, 
                value.var = "residuals")
trainL <- as.matrix(trainL, 1)
somGrid <- somgrid(xdim = 10, ydim = 10, topo= 'hexagonal', toroidal = T)
set.seed(1)
init <- sample(nrow(trainL), somGrid$xdim*somGrid$ydim)
init <- trainL[init, , drop= F]
somL <- supersom(data= trainL, 
                 grid= somGrid, 
                 init= init, 
                 rlen= 100, 
                 keep.data = TRUE,
                 maxNA.fraction = .2)

trainR <- dcast(dat,
                R~L, 
                value.var = "residuals")
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
                 maxNA.fraction = .8)

# Codes clustering
rkcl <- 4
ckcl <- 4
set.seed(7)
clL <- data.table(data.table(somL$data[[1]], keep.rownames = "name"),
                  cl= kmeans(somL$codes[[1]], rkcl)$cluster[somL$unit.classif])
clL[, ord:= median(unlist(.SD[, !"name"]), na.rm = T), cl]
clL[, cl:= .GRP, keyby= ord]
set.seed(5)
clR <- data.table(data.table(somR$data[[1]], keep.rownames = "name"),
                  cl= kmeans(somR$codes[[1]], ckcl)$cluster[somR$unit.classif])
clR[, ord:= median(unlist(.SD[, !"name"]), na.rm = T), cl]
clR[, cl:= .GRP, keyby= ord]

# Heatmap
mat <- trainL[match(clL$name, rownames(trainL)), match(clR$name, colnames(trainL))]
cl <- vl_heatmap(mat,
                 row_clusters = factor(clL$cl, rev(seq(ckcl))),
                 col_clusters = clR$cl,
                 col_clusters_col = rev(grDevices::gray.colors(length(unique(clR$cl)))),
                 cluster_rows = F, 
                 cluster_cols = F,
                 breaks = c(-2,-0.25,0.25,2), 
                 col= c("cornflowerblue", "white", "white", "tomato"), 
                 legend_title = "Obs./Exp. (log2)", 
                 show_rownames = F,
                 show_colnames = F,
                 auto_margins = F,
                 plot= F)
cl$rows[dat, median:= i.median_L, on= "name==L", mult= "first"]
cl$cols[dat, median:= i.median_R, on= "name==R", mult= "first"]

# Check confusion matrix
conf <- merge(cl$rows, cl$cols, by= "name")
tab <- table(conf$cl.x, conf$cl.y)
tab <- matrix(tab, ncol = ncol(tab), dimnames = dimnames(tab))
chi <- chisq.test(tab)
conf <- chi$residuals

# check plot
plot(cl)
# vl_heatmap(x =  matrix(conf, ncol = ncol(conf), dimnames = dimnames(conf)),
#            breaks = c(-3,0,3),
#            cluster_rows = F, 
#            cluster_cols = F, 
#            legend_title = "OR (log2)", 
#            display_numbers = T,
#            display_numbers_matrix = tab)

# saveRDS(cl,
#         "Rdata/vllib002_lm_residuals_SOM.rds")
