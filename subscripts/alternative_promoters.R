setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import ----
dat <- rbindlist(list(dev_highBasal= readRDS("db/FC_tables/vllib028_DESeq2.rds"),
                      dev_lowBasal= readRDS("db/FC_tables/vllib027_DESeq2.rds"),
                      # hk_highBasal= readRDS("db/FC_tables/vllib025_DESeq2.rds"), # Quality is very low
                      hk_lowBasal= readRDS("db/FC_tables/vllib026_DESeq2.rds")),
                 idcol = "class")
dat[, c("class", "cdition"):= tstrsplit(class, "_")]

# Compute add/mult ----
dat[, basalMean:= mean(log2FoldChange[ctlL & ctlR]), .(class, cdition)]
dat[, "Additive model":= log2(2^indL+2^indR-2^basalMean)]
dat[, "Multiplicative model":= log2(2^indL*2^indR/2^basalMean)]
dat[, "Linear model":= (-0.02+1.1*indL+1.1*indR-0.09*indL*indR)]
dat <- dat[, {
  patt <- paste0("^control|^shared|^", class)
  .SD[grepl(patt, L) & grepl(patt, R)]
}, .(class, cdition)]
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]

# Melt ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat,
           id.vars = c("class", "cdition", "Observed", "actPair"),
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl[, Rsq:= {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  1-(((1-rsq)*(.N-1))/(.N-2-1))
}, .(variable, class, cdition)]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)",
                   "Linear model"= "Predicted activity (log2)"), variable]
setorderv(pl,
          c("variable", "class", "cdition"))

# Clip extremes for ploting
clip <- c(0.005, 0.999)

# Plot ----
pdf("pdf/draft/scatterplot_low_high_CPs_models.pdf",
    width = 7.5,
    height = 2.45)
layout(matrix(1:5,
              nrow= 1),
       widths= c(0.5,1,1,1,0.7))
par(font.main= 1,
    las= 1,
    lend= 2,
    cex= 8/12,
    cex.axis= 7/8,
    cex.main= 1,
    mar= c(5,3,2,.5),
    oma= c(0,0,2,0),
    mgp= c(1.5, 0.325, 0),
    tcl= -0.2,
    bty= "n")
pl[, {
  # Performance
  perf <- unique(.SD[, .(variable, Rsq)])
  setorderv(perf, "Rsq")
  perf[, {
    bar <- barplot(Rsq,
                   ylab= "Adj. R2 (goodness of fit)")
    vl_tilt_xaxis(bar, 
                  labels = variable)
    
  }]
  .SD[, {
    # Model
    clip.x <- between(value,
                      quantile(value, clip[1]),
                      quantile(value, clip[2]), incbounds = T)
    clip.y <- between(Observed,
                      quantile(Observed, clip[1]),
                      quantile(Observed, clip[2]), incbounds = T) 
    .SD[clip.x & clip.y, {
      smoothScatter(value,
                    Observed,
                    xlab= xlab,
                    ylab= "Combined activity (log2)",
                    main= variable,
                    xaxt= "n",
                    yaxt= "n",
                    col= adjustcolor(blues9[9], .3))
      
    }]
    axis(1)
    axis(2)
    # Density active pairs
    z <- MASS::kde2d(value[(actPair)],
                     Observed[(actPair)], 
                     h = c(1,1))
    contour(z,
            lwd= 0.25,
            drawlabels= FALSE,
            add= TRUE,
            xpd= NA)
    # R2
    vl_plot_R2(rsquare = Rsq,
               inset= c(-0.1, 0))
    # Abline
    lw <- diff(grconvertX(0:1, "line", "user"))
    lh <- diff(grconvertY(0:1, "line", "user"))
    clip(par("usr")[1]+lw,
         par("usr")[2]-lw,
         par("usr")[3]+lh,
         par("usr")[4]-lh)
    abline(0, 1, lty= "11")
    print(paste0(.GRP, "/", .NGRP))
  }, .(variable, xlab, Rsq)]
  # Boxplot examples
  .c <- c(class, cdition)
  cols <- c(as.character(perf$variable), "Observed")
  dat[(class==.c[1] & cdition==.c[2]) & (actPair), {
    box <- vl_boxplot(.SD,
                      notch = T,
                      # compute_pval= list(c(1,2), c(3,4)),
                      tilt.names= T,
                      ylab= "Activity (log2)")$stats
    lim <- par("usr")[c(3,4)]
    sub <- .SD[apply(.SD, 1, function(x) all(x>lim[1] & x<lim[2]))
               & between(Observed,
                         box[1,4],
                         box[5,4],
                         incbounds = T)]
    sel <- seq(min(sub$Observed),
               max(sub$Observed),
               length.out= 25)
    sub <- lapply(sel, function(x) sub[which.min(abs(Observed-x))])
    sub <- rbindlist(sub)
    sub[, {
      lines(1:4,
            unlist(.BY),
            col= adjustcolor("grey10", .5),
            lwd= 0.5)
      points(1:4,
             unlist(.BY),
             col= adjustcolor("grey10", .5),
             pch= 16,
             cex= 0.4)
    }, (sub)]
    .SD
  }, .SDcols= cols]
  mtext(paste(class, cdition), outer = T)
}, .(class, cdition)]
dev.off()