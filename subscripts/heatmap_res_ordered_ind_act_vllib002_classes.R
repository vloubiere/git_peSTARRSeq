setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
dat[, classL:= tstrsplit(L, "_", keep= 1)]
dat[, classR:= tstrsplit(R, "_", keep= 1)]
dat <- dat[classL==classR]

pdf("pdf/draft/heatmap_residuals_ordered_ind_act_classes.pdf",
    width = 3.85,
    height = 3.25)
layout(matrix(c(2,4,1,3), nrow= 2),
       widths= c(0.1,1),
       heights= c(1,0.1))
par(mar= c(0.2,0.2,0.2, 0.2),
    oma= c(3,3,2,6),
    mgp= c(2, 0.15, 0),
    tcl= -0.1,
    cex.axis= 0.6,
    las= 1,
    xpd= NA)
dat[,{
  # Matrix Inactive candidates
  mat <- dcast(.SD,
               -indL+L~indR+R,
               value.var = "residuals",
               sep = "__")
  mat <- as.matrix(mat[, -1], 1)
  colnames(mat) <- unlist(tstrsplit(colnames(mat), "__", keep= 2))
  while(sum(is.na(mat))>0.05*nrow(mat)*ncol(mat))
  {
    mat <- mat[-which.max(apply(mat, 1, function(x) sum(is.na(x)))),]
    mat <- mat[,-which.max(apply(mat, 2, function(x) sum(is.na(x))))]
  }
  hm <- vl_heatmap(mat, 
                   cluster_rows= F,
                   cluster_cols= F, 
                   breaks = c(-2, 0, 2),
                   show_rownames= F,
                   show_colnames= F, 
                   legend_title= "Residuals (log2)",
                   plot= F,
                   legend.cex = 0.7)
  
  # Heatmap
  plot(hm)
  left <- .SD[rev(hm$rows$name), indL, on= "L", mult = "first"]
  right <- .SD[hm$cols$name, indR, on= "R", mult = "first"]
  title(main= class, line= 1.25)
  title(xlab= "3' activity",
        ylab= "5' activity")
  par(lwd= 0.25)
  # Left
  par(yaxs= "i")
  bar <- barplot(left,
                 col= "grey",
                 space= 0,
                 border= NA,
                 horiz = T,
                 xaxt= "n")
  abline(v= 0, lwd= 0.25, xpd= F)
  axis(side = 3)
  # Right
  par(yaxs= "r",
      xaxs= "i")
  bar <- barplot(right,
                 col= "grey",
                 space= 0,
                 border= NA,
                 yaxt= "n")
  abline(h= 0, lwd= 0.25, xpd= F)
  axis(side = 4)
  plot.new()
  .SD
}, .(class= classL)]
dev.off()