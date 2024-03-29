setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import ----
dat <- rbindlist(list(dev= readRDS("db/FC_tables/vllib006_DESeq2.rds")),
                 idcol = "class")

# Independent activitites
ind <- rbindlist(list(short= readRDS("db/FC_tables/vllib002_DESeq2.rds"),
                      long= readRDS("db/FC_tables/vllib006_DESeq2.rds")),
                 idcol = "spacer")
ind <- rbindlist(list(`5'`= ind[, .(spacer, enh= L, ind= indL)],
                      `3'`= ind[, .(spacer, enh= R, ind= indR)]),
                 idcol= "side")
ind[, side:= factor(side, c("5'", "3'"))]
ind <- dcast(unique(ind), side+enh~spacer, value.var = "ind")
ind <- na.omit(ind)

# Compute add/mult ----
dat[, basalMean:= mean(log2FoldChange[ctlL & ctlR]), class]
dat[, "Additive model":= log2(2^indL+2^indR-2^basalMean)]
dat[, "Multiplicative model":= log2(2^indL*2^indR/2^basalMean)]
dat[, "Linear model":= (-0.02+1.1*indL+1.1*indR-0.09*indL*indR)]
dat <- dat[, {
  patt <- paste0("^control|^shared|^", class)
  .SD[grepl(patt, L) & grepl(patt, R)]
}, class]
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]

# Melt ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat,
           id.vars = c("class", "Observed", "actPair"),
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

# Clip extremes for ploting
clip <- c(0.005, 0.999)

# Plot ----
pdf("pdf/draft/Compare_add_mult_longSpacer_vllib006.pdf",
    width = 10,
    height = 2.45)
layout(matrix(c(1,3,2,4,5,5,6,6,7,7,8,8,9,9),
              nrow= 2),
       widths= c(0.45,0.45,0.45,1,1,1,0.7))
par(bty= "n",
    tcl= -0.1,
    las= 1,
    oma= c(0,0,2,0),
    font.main= 1,
    mgp= c(0.75,0.25,0),
    cex= 0.5,
    cex.axis= 0.5,
    cex.lab= 0.6,
    cex.main= 0.8,
    bty= "n")
# Individual activities
ind[,{
  par(mar= c(4.5,3,1.7,0.2))
  plot(short,
       long,
       frame= F,
       pch= 16,
       col= adjustcolor("grey", 0.4),
       xaxt= "n",
       yaxt= "n",
       xlab= "Act. 300bp spacer (log2)",
       ylab= "Act. 2kb spacer (log2)")
  axis(1, lwd= 0.5)
  axis(2, lwd= 0.5)
  title(main= paste0(side, " cand."),
        xpd= NA,
        line= 0.25)
  PCC <- cor.test(short,
                  long)$estimate
  legend("topleft",
         legend= paste0("PCC= ", round(PCC, 2)),
         bty= "n",
         inset= c(-0.05, 0),
         cex= 0.6)
  abline(0, 1, lty= "11")
  par(mar= c(4.5,4,1.7,1))
  vl_boxplot(short,
             long,
             compute_pval= list(c(1,2)),
             col= "lightgrey",
             ylab= "Individual activity (log2)",
             xaxt= "n",
             yaxt= "n",
             lwd= .5)
  vl_tilt_xaxis(1:2, 
                labels= c("300bp spacer", "2kb spacer"))
  axis(2, lwd= 0.5)
  .SD
}, side]
# Different models
par(lend= 2,
    cex= 8/12,
    cex.axis= 7/8,
    cex.lab= 1,
    cex.main= 1,
    mar= c(5,3,2,.5),
    mgp= c(1.5, 0.325, 0),
    tcl= -0.2)
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
    vl_plot_coeff(value = Rsq,
                  type= "rsq",
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
for(i in 1:5) plot.new()
vl_boxplot(value~variable+side,
           melt(ind, id.vars = "side", measure.vars = c("long", "short")),
           compute_pval= list(c(1,2), c(3,4)),
           col= "lightgrey",
           ylab= "Individual activiy (log2)")
dev.off()