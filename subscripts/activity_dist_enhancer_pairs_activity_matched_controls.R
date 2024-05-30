setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

# Import lib ----
lib <- readRDS("Rdata/vl_library_twist008_112019.rds")
lib <- as.data.table(lib)
lib <- vl_resizeBed(lib, "center", 0, 0)

# Import dat and compute left right distance ----
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
dat[lib, c("seqL", "startL"):= .(i.seqnames, i.start), on="L==ID_vl"]
dat[lib, c("seqR", "startR"):= .(i.seqnames, i.start), on="R==ID_vl"]
dat[, dist:= ifelse(seqL!=seqR, Inf, abs(startL-startR))]
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Find best activity-matched controls based on indL and indR----
dat$ctlL <- dat$ctlR <- NULL
enr <- dat[dist<=20e3]
ctls <- dat[dist>20e3]
# To limit search space, max diff between enr and ctls
breaks <- 0.1
enr[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
ctls[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
# While some pairs have no control
enr[, ctlL:= as.character(NA)]
tot <- nrow(enr)
while(nrow(ctls) && any(is.na(enr$ctlL)))
{
  # For each enr pair, define the ctl pair with smallest euclidean distance
  .c <- ctls[enr[is.na(ctlL)], {
    dist <- sqrt((x.indL-i.indL)^2+(x.indR-i.indR)^2)
    idx <- which.min(dist)
    .(L= i.L,
      R= i.R,
      ctlL= x.L[idx],
      ctlIndL= x.indL[idx],
      ctlR= x.R[idx],
      ctlIndR= x.indR[idx],
      dist= dist[idx])
  }, .EACHI, on= c("breakL", "breakR"), nomatch= NULL]
  # for each (potentially duplicated) control pair, select closest enr pair
  setorderv(.c, "dist")
  .c <- .c[, .SD[1], .(ctlL, ctlR)]
  # Add control pairs to enr
  enr[.c, c("ctlL", "ctlIndL", "ctlR", "ctlIndR"):= .(i.ctlL, i.ctlIndL, i.ctlR, i.ctlIndR), on= c("L", "R")]
  # Remove control pairs from potential controls
  ctls[.c, used:= T, on= c("L==ctlL", "R==ctlR")]
  ctls <- ctls[is.na(used), !"used"]
  # Increase max dist when no more control pairs can be found
  if(!nrow(.c))
  {
    breaks <- breaks+0.5
    enr[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
    ctls[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
  }
  print(paste0(sum(!is.na(enr$ctlL)), "/", tot))
}
enr[dat, ctlLog2FoldChange:= i.log2FoldChange, on= c("ctlL==L", "ctlR==R")]

# Plot ----
Cc <- c("lightgrey", "rosybrown1")
pdf("pdf/draft/closeby_enhancers_vs_distant_controls.pdf", 
    width = 3, 
    height = 3)
par(mai= c(0.75,0.85,1.2,1.5), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75)
xpos <- seq(1, 5, length.out= 3)
vl_boxplot(enr[, .(ctlIndL, indL, ctlIndR, indR, ctlLog2FoldChange, log2FoldChange)],
           compute.pval= list(c(1,2), c(3,4), c(5,6)),
           notch= T,
           xaxt= "n",
           col= Cc,
           ylab= "Activity (log2)",
           at= rep(xpos, each= 2)+c(-0.35, 0.35),
           lwd= .5)
axis(1, at= xpos, labels = c("5'", "3'", "Pair"))
legend(par("usr")[2],
       par("usr")[4],
       fill= Cc,
       legend= c(">20kb in situ",
                 "<20kb in situ"),
       bty= "n",
       xpd= NA,
       cex= 0.7)
dev.off()