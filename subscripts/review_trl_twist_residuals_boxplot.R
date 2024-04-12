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

pl <- list("All pairs"= dat,
           "Twist/Twist"= dat[Twist__L>0 & Twist__R>0],
           "GATA/GATA"= dat[GATA__L>0  & GATA__R>0])
.n <- c("All pairs", "Twist/Twist", "GATA/GATA")
names(pl)
pl <- rbindlist(pl, idcol = T)
pl[, .id:= factor(.id, .n)]
  
# Plot ----
Cc <- c("lightgrey", "limegreen", "tomato")
pdf("pdf/draft/review_Trl_Twist_pairs_residuals_boxplot.pdf", 5, 3)
vl_par(mfrow= c(1, 2),
       mgp= c(1, .35, 0))
box <- vl_boxplot(residuals~.id, 
                  pl,
                  names= levels(pl$.id), 
                  compute.pval= list(c(1,2),
                                     c(1,3)),
                  tilt.names= T,
                  ylab= "Residuals (log2)",
                  col= Cc)
abline(h= box$stats[3,c(1, 2)],
       lty= "11",
       col= Cc[1:2])
var <- round(box$stats[3,c(1, 2)], 1)
var <- formatC(var, format = "f", digits = 1)
text(par("usr")[2]+strwidth(var, cex= .7)/2,
     box$stats[3,c(1, 2)],
     var,
     pos= c(1,3),
     cex= .7,
     xpd= NA,
     offset= 0,
     col= Cc[c(1,2)])
box <- vl_boxplot(log2FoldChange~.id, 
                  pl,
                  names= levels(pl$.id), 
                  compute.pval= list(c(1,2),
                                     c(1,3)),
                  tilt.names= T,
                  ylab= "Activity (log2)",
                  col= Cc)
abline(h= box$stats[3,c(1, 3)],
       lty= "11",
       col= Cc[c(1,3)])
var <- round(box$stats[3,c(1, 3)], 1)
var <- formatC(var, format = "f", digits = 1)
text(par("usr")[2]+strwidth(var, cex= .7)/2,
     box$stats[3,c(1, 3)],
     var,
     pos= c(1,3),
     cex= .7,
     xpd= NA,
     offset= 0,
     col= Cc[c(1,3)])
dev.off()