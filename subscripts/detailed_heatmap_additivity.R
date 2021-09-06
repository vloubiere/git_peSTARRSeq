setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
if("plot_group" %in% names(dat))
  dat$plot_group <- NULL
# dat[cdition=="vllib002" & ((group_L=="OSC" & group_R=="OSC" & median_L>1 & median_R>1)), plot_group:= "vllib002: OSC x OSC"]
# dat[cdition=="vllib002" & ((group_L=="OSC" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: OSC x dev"]
# dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x OSC"]
# dat[cdition=="vllib002" & ((group_L=="hk" & group_R=="hk" & median_L>1 & median_R>1)), plot_group:= "vllib002: hk x hk"]
# dat[cdition=="vllib002" & ((group_L=="hk" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: hk x dev"]
# dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="hk" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x hk"]
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x dev"]

# dat[cdition=="vllib015" & group_L=="hk" & group_R=="hk", plot_group:= "vllib015: hk x hk"]
# dat[cdition=="vllib015" & group_L=="hk" & group_R=="dev", plot_group:= "vllib015: hk x dev"]
# dat[cdition=="vllib015" & group_L=="dev" & group_R=="hk", plot_group:= "vllib015: dev x hk"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="dev", plot_group:= "vllib015: dev x dev"]
# dat[cdition=="vllib015" & group_L=="DHS" & group_R=="DHS", plot_group:= "vllib015: DHS x DHS"]
# dat[cdition=="vllib015" & group_L=="dev" & group_R=="DHS", plot_group:= "vllib015: dev x DHS"]
# dat[cdition=="vllib015" & group_L=="DHS" & group_R=="dev", plot_group:= "vllib015: DHS x dev"]
# dat[cdition=="vllib015" & group_L=="Silencer" & group_R=="Silencer", plot_group:= "vllib015: Sil. x Sil."]
# dat[cdition=="vllib015" & group_L=="Silencer" & group_R=="dev", plot_group:= "vllib015: Sil. x dev"]
# dat[cdition=="vllib015" & group_L=="dev" & group_R=="Silencer", plot_group:= "vllib015: dev x Sil."]
# dat[cdition=="vllib015" & group_L=="CP" & group_R=="CP", plot_group:= "vllib015: CP x CP"]
# dat[cdition=="vllib015" & group_L=="CP" & group_R=="dev", plot_group:= "vllib015: CP x dev"]
# dat[cdition=="vllib015" & group_L=="dev" & group_R=="CP", plot_group:= "vllib015: dev x CP"]
# dat[cdition=="vllib015" & group_L=="SUHW_peak" & group_R=="SUHW_peak", plot_group:= "vllib015: Ins. x Ins."]
# dat[cdition=="vllib015" & group_L=="SUHW_peak" & group_R=="dev", plot_group:= "vllib015: Ins. x dev"]
# dat[cdition=="vllib015" & group_L=="dev" & group_R=="SUHW_peak", plot_group:= "vllib015: dev x Ins."]
# dat[cdition=="vllib015" & group_L=="SUHW_peak" & group_R=="SUHW_peak", plot_group:= "vllib015: Ins. x Ins."]

dat[cdition=="vllib016" & group_L=="hk" & group_R=="hk", plot_group:= "vllib016: hk x hk"]
# dat[cdition=="vllib016" & group_L=="hk" & group_R=="dev", plot_group:= "vllib016: hk x dev"]
# dat[cdition=="vllib016" & group_L=="dev" & group_R=="hk", plot_group:= "vllib016: dev hk"]
# dat[cdition=="vllib016" & group_L=="dev" & group_R=="dev", plot_group:= "vllib016: dev x dev"]
# dat[cdition=="vllib016" & group_L=="DHS" & group_R=="DHS", plot_group:= "vllib016: DHS x DHS"]
# dat[cdition=="vllib016" & group_L=="DHS" & group_R=="hk", plot_group:= "vllib016: DHS x hk"]
# dat[cdition=="vllib016" & group_L=="hk" & group_R=="DHS", plot_group:= "vllib016: hk x DHS"]
# dat[cdition=="vllib016" & group_L=="Silencer" & group_R=="Silencer", plot_group:= "vllib016: Sil. x Sil."]
# dat[cdition=="vllib016" & group_L=="Silencer" & group_R=="hk", plot_group:= "vllib016: Sil. x hk"]
# dat[cdition=="vllib016" & group_L=="hk" & group_R=="Silencer", plot_group:= "vllib016: hk x Sil."]
# dat[cdition=="vllib016" & group_L=="SUHW_peak" & group_R=="SUHW_peak", plot_group:= "vllib016: Ins. x Ins."]
# dat[cdition=="vllib016" & group_L=="SUHW_peak" & group_R=="hk", plot_group:= "vllib016: Ins. x hk"]
# dat[cdition=="vllib016" & group_L=="hk" & group_R=="SUHW_peak", plot_group:= "vllib016: hk x Ins."]
# dat[cdition=="vllib016" & group_L=="CP" & group_R=="CP", plot_group:= "vllib016: CP x CP"]
# dat[cdition=="vllib016" & group_L=="CP" & group_R=="hk", plot_group:= "vllib016: CP x hk"]
# dat[cdition=="vllib016" & group_L=="hk" & group_R=="CP", plot_group:= "vllib016: hk x CP"]


pdf("pdf/all_pairs_activity_additivity.pdf", width = 8.2)
layout(matrix(c(2,1,4,3), ncol= 2), 
       widths = c(1, 0.15),
       heights = c(0.15, 1))

dat[!is.na(plot_group) & L!=R, {
  # Use expected additive values to order the heatmap
  .mat <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "diff")
  .mat <- as.matrix(.mat, 1)
  row_ord <- order(-.SD[rownames(.mat), median_L[1], L, on= "L"]$V1)
  col_ord <- order(.SD[colnames(.mat), median_R[1], R, on= "R"]$V1)
  .mat <- .mat[row_ord, ]
  .mat <- .mat[, col_ord]
  
  # Margins
  mBot= 0.5
  mLeft= 0.5
  mTop= 0
  mRight= 0.98
  par(mai= c(mBot, mLeft, mTop, mRight), 
      mgp = c(3, 0.4, 0))
  
  # Heatmap
  .breaks <- max(abs(quantile(.mat, c(0.05, 0.975), na.rm= T)))
  res <- vl_heatmap(mat = .mat, 
                    # main= plot_group, 
                    cluster_rows = F,
                    cluster_cols = F,
                    show_rownames = F,
                    show_colnames = F,
                    breaks = c(-.breaks, -0.5, 0.5, .breaks),
                    col= c("cornflowerblue", "white", "white", "tomato"),
                    legend_title = "Additive Score (log2)", 
                    auto_margins = F)
  rect(0, 0, 1, 1, lwd= 0.5)
  mtext("Right Regulatory Element", 
        1, 
        cex= 2)
  mtext("Left Regulatory Element", 
        2, 
        cex= 2)
  
  # Median activities
  for(side in c("right", "left"))
  {
    if(side=="right")
    {
      .c <- .SD[colnames(.mat), .(value= median_R[1], 
                                  col= col_R[1], 
                                  y= .GRP), .(idx= R), on= "R"]
      .f <- x~y
      axPos <- 2
      axLine <- 0.6
      axLas <- 1
      .liney <- vl_scale01(colMeans(.mat, na.rm= T))
      .linex <- seq(0, 1, length.out = length(.liney))
    }
    if(side=="left")
    {
      .c <- .SD[rev(rownames(.mat)), .(value= median_L[1], 
                                  col= col_L[1], 
                                  y= .GRP), .(idx=L), on= "L"]
      .f <- y~x
      axPos <- 1
      axLine <- 1.5
      axLas <- 0
      .linex <- vl_scale01(rowMeans(.mat, na.rm= T))
      .liney <- seq(1, 0, length.out = length(.linex))
    }
      
    # Compute im
    .c <- .c[, .(x= 1:100,
                 Cc= "white"), (.c)]
    .c[, max:= round(value/max(value)*100), ]
    .c[x<=max, Cc:= col]
    im <- dcast(.c,
                .f, 
                value.var = "Cc")
    im <- as.matrix(im, 1)
    par(mai= c(ifelse(side=="right", 0, mBot),
               ifelse(side=="right", mLeft, 0),
               ifelse(side=="right", 0, mTop),
               ifelse(side=="right", mRight, 0)))
    # Plot image
    plot.new()
    rasterImage(im, 
                interpolate = F,
                xleft = 0, 
                xright = 1,
                ytop = 0, 
                ybottom = 1)
    ticks <- axisTicks(c(0, max(.c$value)), log= F)
    axis(axPos,
         labels= ticks,
         at= ticks/max(.c$value),
         xpd= T,
         line = -1, 
         tcl= -0.25,
         cex= 0.6,
         las= axLas)
    mtext("Individual\nActivity (log2)", 
          axPos, 
          cex= 0.8,
          line = axLine)
    # Add Marginals
    lines(.linex, 
          .liney)
  }
  plot.new()
  
  print("DONE")
}, plot_group]
dev.off()


