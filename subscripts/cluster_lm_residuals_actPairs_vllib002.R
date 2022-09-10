setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(kohonen)

# Import data
dat <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh."]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
dat[, predicted:= predict(model, newdata = .SD)]
dat[, residuals:= log2FoldChange-predicted]

###############################################
# CLUSTERING
###############################################
if(!file.exists("Rdata/clustering_lm_residuals_vllib002.rds"))
{
  # SOM clustering
  trainL <- dcast(dat,
                  L~R, 
                  value.var = "residuals")
  trainL <- as.matrix(trainL, 1)
  # clip <- quantile(trainL, c(0.001, 0.999), na.rm= T)
  # trainL[trainL<clip[1]] <- clip[1]
  # trainL[trainL>clip[2]] <- clip[2]
  trainL <- trainL[order(-dat[rownames(trainL), median_L, on= "L", mult= "first"]),]
  somGrid <- somgrid(xdim = 10, ydim = 10, topo= 'hexagonal', toroidal = T)
  set.seed(1)
  init <- sample(nrow(trainL), somGrid$xdim*somGrid$ydim)
  init <- trainL[init, , drop= F]
  somL <- supersom(data= trainL, 
                   grid= somGrid, 
                   init= init, 
                   rlen= 100, 
                   keep.data = TRUE,
                   maxNA.fraction = .3)
  
  trainR <- dcast(dat,
                  R~L, 
                  value.var = "residuals")
  trainR <- as.matrix(trainR, 1)
  # clip <- quantile(trainR, c(0.001, 0.999), na.rm= T)
  # trainL[trainR<clip[1]] <- clip[1]
  # trainL[trainR>clip[2]] <- clip[2]
  trainR <- trainR[order(dat[rownames(trainR), median_R, on= "R", mult= "first"]),]
  set.seed(1)
  init <- sample(nrow(trainR), somGrid$xdim*somGrid$ydim)
  init <- trainR[init, , drop= F]
  somR <- supersom(data= trainR, 
                   grid= somGrid, 
                   init= init, 
                   rlen= 100, 
                   keep.data = TRUE,
                   maxNA.fraction = .3)
  
  # Codes clustering
  par(mfrow= c(2,1))
  set.seed(1)
  wss <- sapply(1:15, function(k){kmeans(somL$codes[[1]], k, nstart=50, iter.max = 15)$tot.withinss})
  plot(wss, xlab= "Left Number of clusters", ylab= "wss")
  lines(wss)
  set.seed(1)
  wss <- sapply(1:15, function(k){kmeans(somR$codes[[1]], k, nstart=50, iter.max = 15)$tot.withinss})
  plot(wss, xlab= "Left Number of clusters", ylab= "wss")
  lines(wss)
  # set.seed(2)
  # somL$unit.classif.kcl <- kmeans(somL$codes[[1]], 4)$cluster[somL$unit.classif]
  somL$unit.classif.kcl <- cutree(hclust(dist(somL$codes[[1]])), 4)[somL$unit.classif]
  # set.seed(6)
  # somR$unit.classif.kcl <- kmeans(somR$codes[[1]], )$cluster[somR$unit.classif]
  somR$unit.classif.kcl <- cutree(hclust(dist(somR$codes[[1]])), 4)[somR$unit.classif]
  
  # Heatmap
  mat <- somL$data[[1]]
  mat <- mat[, match(rownames(somR$data[[1]]), colnames(mat))]
  cl <- vl_heatmap(mat,
                   row_clusters = factor(somL$unit.classif.kcl, 
                                         labels = c("High syn.", "Medium syn.", "Weak syn.", "No syn.")),
                   col_clusters = factor(somR$unit.classif.kcl,
                                         labels = c("No syn.", "Weak syn.", "Medium syn.", "High syn."),
                                         levels = c(2,3,1,4)),
                   col_clusters_col = grey.colors(5)[-1],
                   row_clusters_col = rev(grey.colors(5)[-1]),
                   cluster_rows = F, 
                   cluster_cols = F,
                   breaks = c(-1.5,-0.25, 0.25, 1.5), 
                   col= c("cornflowerblue", "white", "white", "tomato"), 
                   legend_title = "Obs./Exp. (log2)", 
                   show_rownames = F,
                   show_colnames = F,
                   plot= F)
  cl$rows[dat, c("median", "class_act"):= .(i.median_L, i.class_act_L), on= "name==L", mult= "first"]
  cl$cols[dat, c("median", "class_act"):= .(i.median_R, i.class_act_R), on= "name==R", mult= "first"]
  
  # SAVE
  saveRDS(cl, 
          "Rdata/clustering_lm_residuals_vllib002.rds")
}

###############################################
# PLOT Heatmap
###############################################
cl <- readRDS("Rdata/clustering_lm_residuals_vllib002.rds")

pdf("pdf/draft/heatmap_som_clustering_residuals_vllib002.pdf")
layout(matrix(c(2,4,1,3), ncol= 2),
       widths = c(0.6,10),
       heights = c(10,0.6))
par(mar= c(0.2,0.2,3,5),
    tcl= -0.2,
    las= 1)
plot(cl)
par(mar= c(0.2,0.4,3,0.5),
    yaxs= "i",
    cex.axis= 0.5,
    mgp= c(2,0.1,0.2),
    tcl= -0.15)
cl$rows[order(y), {
  barplot(median, 
          horiz= T,
          border= NA, 
          col= "grey30",
          space= 0, 
          xaxt= "n",
          xlim= c(0, 8.5))
  abline(v= 0, lwd= 0.25)
  axis(3, at= c(0, 8), labels = c(0, 8))
  text(4.5,
       par("usr")[4],
       "5' act\n(log2)",
       xpd= T, 
       pos= 3, 
       offset= 1.15,
       cex= 0.6)
}]
par(mar= c(0.4,0.2,0.5,5),
    yaxs= "r",
    xaxs= "i",
    mgp= c(2,0.2,0.2),
    tcl= -0.1)
cl$cols[order(x), {
  barplot(median,
          border= NA, 
          col= "grey30",
          space= 0, 
          yaxt= "n",
          ylim= c(0, 8.5))
  abline(h= 0, lwd= 0.25)
  axis(4, at= c(0, 8), labels = c(0, 8))
  text(par("usr")[2],
       3,
       "3' act\n(log2)",
       xpd= T,
       pos= 4,
       offset= 0.8,
       cex= 0.6)
}]
dev.off()