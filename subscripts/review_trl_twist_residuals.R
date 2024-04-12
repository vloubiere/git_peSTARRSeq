setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")

# Import motifs ----
motID <- readRDS("db/motif_counts/lib8_motifs_IDs.rds")
motID <- motID[, lapply(.SD, as.character)]
motID[, name:= switch(cluster,
                      "Ebox/CATATG/twi"= "Twist",
                      "Trl/1"= "Trl", 
                      "AP1/1"= "AP-1", 
                      "GATA/1"= "GATA"), cluster]
motID <- motID[!is.na(name)]
mot <- readRDS("db/motif_counts/twist008_motif_counts.rds")
mot <- mot[, c("ID", motID$motif_ID), with= F]
setnames(mot, names(mot)[-1], motID$name)
dat <- merge(dat, mot, by.x= "L", by.y= "ID")
dat <- merge(dat, mot, by.x= "R", by.y= "ID", suffixes= c("__L", "__R"))

# Merge ----
left <- dat[, .(res= mean(residuals), act= mean(log2FoldChange)), .(enh= L, Twist= Twist__L, Trl= Trl__L)]
right <- dat[, .(res= mean(residuals), act= mean(log2FoldChange)), .(enh= R)]
.m <- merge(left, right, by= "enh", suffixes= c("_L", "_R"))
.m[, Cc:= ifelse(Twist>0, "limegreen", "lightgrey")]
.m[, col:= ifelse(Trl>0, "black", NA)]
.m <- .m[order(Cc!="lightgrey")]

# Ploting function ----
pl <- function(x= "res_R",
               y= "res_L",
               xlab= "Mean residuals 3'",
               ylab= "Mean residuals 5'",
               ab= TRUE)
{
  .m[, {
    plot(get(x),
         get(y),
         pch= 21,
         bg= adjustcolor(Cc, .5),
         col= col,
         cex= .5,
         xlab= xlab,
         ylab= ylab,
         lwd= .35,
         xaxt= "n")
    axis(1, padj= -1.25)
  }]
  if(ab)
  {
    abline(h= 0, lty= "11")
    abline(v= 0, lty= "11")
  }
  dens <- list(all= .m[Trl==0 & Twist==0],
               Twist= .m[Twist>0],
               Trl= .m[Trl>0])
  dens <- rbindlist(dens, idcol = T)
  dens[, Cc:= c("lightgrey", "limegreen", "black")[.GRP], .id]
  yadj <- 
  dens[, {
    .d <- density(get(x), from= par("usr")[1], to= par("usr")[2])
    .d$y <- .d$y/max(.d$y)*strheight("M")*1.5+par("usr")[4]
    lines(.d, xpd= T, col= Cc[1], lwd= .75)
    segments(.d$x[which.max(.d$y)],
             par("usr")[4],
             .d$x[which.max(.d$y)],
             max(.d$y),
             xpd= NA,
             col= Cc[1],
             lwd= .5,
             lty= "11")
    .d <- density(get(y), from= par("usr")[3], to= par("usr")[4])
    .d$y <- .d$y/max(.d$y)*strwidth("M")*1.5+par("usr")[2]
    lines(.d$y, .d$x, xpd= T, col= Cc[1], lwd= .75)
    segments(par("usr")[2],
             .d$x[which.max(.d$y)],
             max(.d$y),
             .d$x[which.max(.d$y)],
             xpd= NA,
             col= Cc[1],
             lwd= .5,
             lty= "11")
  }, .(.id, Cc)]
  legend("topleft",
         bty= "n",
         pt.bg= c("lightgrey", "white", "limegreen"),
         pt.lwd= .35,
         legend= c("No Trl/Twist", "Trl", "Twist"),
         pch= 21,
         col= c(NA, "black", NA),
         cex= .5)
}

# Plot ----
pdf("pdf/draft/review_Trl_Twist_pairs_residuals.pdf", 6, 3)
vl_par(mfrow= c(1,2))
pl()
pl("act_R", "act_L", "Mean activity 3'", "Mean activity 5'", ab= F)
dev.off()