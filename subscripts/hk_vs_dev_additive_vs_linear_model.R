setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import ----
dat <- rbindlist(list(dev= readRDS("db/FC_tables/vllib015_DESeq2.rds"), # DSCP core promoter
                      hk= readRDS("db/FC_tables/vllib016_DESeq2.rds")), # RpS12 core promoter
                 idcol = "class")
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Compute additive / multiplicative / linear models ----
dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]
dat[, `Linear model`:= predict(readRDS("db/linear_models/lm_vllib002.rds"), .SD)]
dat <- dat[, { # Select dev/hk enhancer-enhancer pairs
  patt <- paste0("^control|^", class)
  .SD[grepl(patt, L) & grepl(patt, R)]
}, class]
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]

# Melt data ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat,
           id.vars = c("class", "L", "R", "Observed", "actPair"),
           measure.vars = c("Additive model", "Multiplicative model", "Linear model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl[, Rsq:= {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  1-(((1-rsq)*(.N-1))/(.N-2-1))
}, .(variable, class)]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)",
                   "Linear model"= "Predicted activity (log2)"), variable]
setorderv(pl,
          c("variable", "class"))

# Check best predictor
t1 <- pl[(actPair) & variable!="Linear model" & class=="hk", .(variable[which.min(abs(residuals))]), .(L, R)]$V1
t1 <- table(droplevels(t1))

# Clip extremes for ploting
clip <- c(0.005, 0.999)

# Plot ----
pdf("pdf/draft/scatterplot_dev_hk_models.pdf",
    width = 15,
    height = 3)
par(mgp= c(0.75, 0.25, 0),
    mfrow= c(1,5),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
pl[, {
  par(mai= rep(.9, 4), 
      pty= "s")
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
                    main= paste(variable, class),
                    xaxt= "n",
                    yaxt= "n",
                    col= adjustcolor(blues9[9], .3))
      
    }]
    axis(1, padj= -1.25)
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
    leg.x <- par("usr")[1]+strwidth("M")*0.75
    leg.y <- par("usr")[4]-strheight("M")*.75
    points(x = rep(leg.x, 2), # Density legend
           y = rep(leg.y, 2),
           lwd= 0.5,
           cex= c(0.6, 1.25))
    text(x = leg.x,
         y = leg.y,
         pos= 4,
         labels= "enh./enh. pairs",
         cex= 7/12)
    # R2
    vl_plot_coeff(value = Rsq,
                  inset= c(-0.075, 0.06),
                  cex= 7/12)
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
  # Performances
  par(mai= c(.9,1.25,.9,1.25),
      pty= "m",
      mgp= c(1, .25, 0))
  perf <- unique(.SD[, .(variable, Rsq)])
  setorderv(perf, "Rsq")
  perf[, {
    bar <- barplot(Rsq,
                   ylab= "Adj. R2 (goodness of fit)")
    vl_tilt_xaxis(bar, 
                  labels = variable)
    
  }]
  # Boxplot examples
  par(mai= c(.9,1.25,.9,1.1),
      mgp= c(.75, .25, 0))
  .c <- class
  cols <- c(as.character(perf$variable), "Observed")
  dat[(class==.c) & (actPair), {
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
  mtext(class, outer = T)
}, class]
plot.new()
par(mai= rep(1, 4),
    pty= "s")
pie(t1,
    labels = paste0(names(t1), " (", round(t1/sum(t1)*100), "%) n= ", formatC(t1, big.mark = ",")),
    cex= 8/12,
    lwd= .5,
    xpd= NA)
dev.off()