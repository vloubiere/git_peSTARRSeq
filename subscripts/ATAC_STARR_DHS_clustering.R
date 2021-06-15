dat <- fread("Rdata/STARR_ATAC_DHS_CHIP_quantif.txt")

#---------------------------#
# Clustering
#---------------------------#
mat <- as.matrix(dat[, -c(1:5)])
rownames(mat) <- paste0(dat$seqnames, ":", dat$start, "-", dat$end)

require(kohonen)
grid <- somgrid(xdim= 4, 
                ydim = 4, 
                topo = "hexagonal", 
                toroidal = T)
set.seed(1)
# init <- mat[sample(nrow(mat), grid$xdim*grid$ydim),]
# model <- supersom(scale(mat), 
#                   grid,
#                   init = init)

model <- supersom(list(mat[, 1:4],
                       mat[, 5:8]),
                  grid, user.weights = c(1, 5))

saveRDS(model, "Rdata/STARR_ATAC_DHS_CHIP_SOM_clustering.rds")

#---------------------------#
# Plots
#---------------------------#
# plot(model, shape= "straight")

#1 ------ HEATMAP ------ ######
pdf("pdf/clustering_ATAC_DHS_STARR.pdf", width = 4.5)
par(mar= c(8,5,2,10), xaxs= "i", yaxs= "i")
cl <- vl_heatmap(mat[order(model$unit.classif),], 
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = F,
                 cutree_rows = 5, 
                 col = c("cornflowerblue", "white", "white", "white", "tomato"),
                 breaks = c(-5,-0.5,0,0.5,5), 
                 legend_title = "Enrichment (log2)")

#Add row cluster to clustering object
cl[data.table(row= rownames(model$data[[1]]), cl= model$unit.classif), rcl:= i.cl, on= "row"]
at <- cl[, .(V1= max(y)/nrow(mat), N= length(unique(row))), rcl]
abline(h= at$V1[-1])
at[, V2:= c(V1[-1], 0)]
axis(2, 
     at= rowMeans(at[, .(V1, V2)]), 
     labels = paste0(at$rcl, " (", at$N, ")"), 
     lwd= 0,
     las= 1, 
     cex.axis= 0.7)
dev.off()

#1 ------ SAVE CKUSTERING OBJECT ------ ######
fwrite(cl, 
       "Rdata/STARR_ATAC_DHS_cluster_final.txt", 
       col.names = T, 
       row.names = F, 
       sep= "\t", 
       quote= F)



