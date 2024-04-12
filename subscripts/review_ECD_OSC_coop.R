setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- list("S2"= readRDS("db/FC_tables/DSCP_ECD_WT_DESeq2.rds")[(grepl("^dev", L) & grepl("^dev", R))],
            "S2+ECD"= readRDS("db/FC_tables/DSCP_ECD_WT_DESeq2.rds")[(grepl("^ecd", L) & grepl("^ecd", R))],
            "OSC"= readRDS("db/FC_tables/DSCP_OSC_WT_DESeq2.rds")[(grepl("^dev", L) & grepl("^dev", R))])
dat <- rbindlist(dat, idcol = "cdition")

# Classes ----
dat[, cdition:= factor(cdition, c("S2", "S2+ECD", "OSC"))]
dat[, classL:= tstrsplit(L, "_", keep= 1), L]
dat[, classR:= tstrsplit(R, "_", keep= 1), R]
dat[, actPair:= actL!="Inactive" & actR!="Inactive"]
dat <- dat[(actPair)]

# Compute expected scores ----
lm <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
dat[, expected:= predict(lm, dat)]
dat[, residuals:= log2FoldChange-(indL+indR)]
dat[, residuals:= log2FoldChange-log2(2^indL+2^indR-1)]

Cc <- c("white", "cadetblue2", "tan1")
par(mfrow= c(3,1))
dat[, {
  plot(NA,
       xlim= c(-5, 5),
       ylim= c(0, 0.6),
       xlab= "Residuals (log2)",
       ylab= "Density")
  hist(residuals,
       freq = F, 
       add= T,
       border= adjustcolor("black", .3),
       col= adjustcolor(Cc[.GRP], .3))
  print("")
}, cdition]

hist(dat$residuals)

dat[, `Additive model`:= log2(2^indL+2^indR-1)]
dat[, `Multiplicative model`:= indL+indR]


dat[cdition=="OSC", vl_rasterScatterplot(`Additive model`, log2FoldChange, col= ifelse(grepl("^dev", L), "green", "red"))]

# Melt data for plotting ----
setnames(dat, "log2FoldChange", "Observed")
pl <- melt(dat, 
           id.vars = c("cdition", "L", "R", "Observed", "actPair"), 
           measure.vars = c("Additive model", "Multiplicative model"))
pl <- na.omit(pl)
pl[, residuals:= Observed-value]
pl[, Rsq:= {
  rsq <- vl_model_eval(Observed, value)$Rsquare
  1-(((1-rsq)*(.N-1))/(.N-2-1))
}, .(variable, cdition)]
pl[, xlab:= switch(as.character(variable), 
                   "Additive model"= "5' + 3' activities (log2)",
                   "Multiplicative model"= "5' x 3' activities (log2)"), variable]

# Check best predictor ----
t1 <- pl[(actPair), .(variable[which.min(abs(residuals))]), .(L, R)]$V1
t1 <- table(droplevels(t1))
t2 <- pl[(actPair), .(variable[which.min(abs(residuals))]), .(L, R)]$V1
t2 <- table(droplevels(t2))

# Clip extremes for ploting ----
clip <- c(0.005, 0.999)

pdf("pdf/draft/Compare_add_mult_ECD_OSC.pdf",
    width = 15,
    height = 3)
layout(matrix(1:5,
              nrow= 1))
par(mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
pl[, {
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
  title(sub = cdition)
  axis(1, padj = -1.25)
  axis(2)
  # Density active pairs
  z <- MASS::kde2d(value[(actPair)],
                   Observed[(actPair)], 
                   n = 50)
  contour(z,
          lwd= 0.25,
          drawlabels= FALSE,
          nlevels= 6,
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
}, .(variable, xlab, Rsq, cdition)]
dev.off()