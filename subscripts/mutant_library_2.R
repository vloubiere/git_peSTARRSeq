setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_mutant_library_FC_DESeq2.rds")

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("groupL", "mutL", "L"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "R"):= tstrsplit(R, "_", keep= c(1,2,4))]
dat[, L:= paste0(groupL, L)]
dat[, R:= paste0(groupR, R)]
dat <- dat[!(grepl("^control", L) & mutL=="WT") & !(grepl("^control", R) & mutR=="WT")]

# Compute expected additive and fitted multiplicative ----
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
dat[, `linear model`:= predict(model, .SD)]
dat[, `additive model`:= log2(2^indL+2^indR-1)]
dat[, residuals:= log2FoldChange-`additive model`]
dat <- dat[actL!="Inactive" & actR!="Inactive"]

# Fitted curves ----
mult <- data.table(indL= seq(0, 10, .1))
mult[, indR:= indL]
mult[, add:= log2(2^indL+2^indR-1)]
mult[, `mult`:= predict(model, .SD)]

# Merge WT and mutant variants ----
pl <- merge(dat[mutL=="WT" & mutR=="WT", .(L, R, indL, indR, log2FoldChange, `linear model`, `additive model`, residuals)],
            dat[mutL!="WT" & mutR!="WT", .(L, R, mutL, mutR, indL, indR, log2FoldChange, `linear model`, `additive model`, residuals)],
            by= c("L", "R"),
            suffixes= c("_wt", "_mut"))
# Retrieve combinations of interest
pl[, cdition:= fcase(grepl("add.*Trl", mutL) & grepl("add.*Trl", mutR), "Added Trl motifs",
                     grepl("add.*Twist", mutL) & grepl("add.*Twist", mutR), "Added Twist motifs",
                     grepl("add.*Dref", mutL) & grepl("add.*Dref", mutR), "Added Dref motifs",
                     grepl("mut.*Trl", mutL) & grepl("mut.*Trl", mutR), "Mutated Trl motifs",
                     grepl("mut.*Twist", mutL) & grepl("mut.*Twist", mutR), "Mutated Twist motifs")]
pl <- pl[!is.na(cdition)]

# Plotting function adds legend to residuals boxplot ----
leg <- function(cdition, n)
{
  var <- fcase(grepl("Mutated", cdition), "Mutant",
               grepl("Added", cdition), "Added motif")
  var <- paste0(var, " (n= ", n, ")")
  legend(par("usr")[1],
         par("usr")[4],
         col= rev(Cc),
         pch= 16,
         legend= c(var, "WT"),
         bty= "n",
         xpd= NA,
         cex= 6/12,
         x.intersp= 0.5)
}

# Plot ----
pdf("pdf/draft/mutant_library_per_cdition.pdf", 
    width = 9,
    height= 15)
par(mfrow= c(5,3),
    mgp= c(.75, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    font.main= 1,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75)
Cc <- c("grey40", "pink1")
pl[, {
  par(pty="s",
      mai= rep(.9 ,4))
  # Observed vs expected
  x <- c(`additive model_wt`, `additive model_mut`)
  y <- c(`log2FoldChange_wt`, `log2FoldChange_mut`)
  col <- rep(Cc, each= .N)
  set.seed(1)
  rdm <- sample(.N*2, .N*2)
  plot(x[rdm],
       y[rdm], 
       main= cdition,
       cex= .5,
       pch= 16,
       col= col[rdm],
       xlab= "Predicted additive (log2)",
       ylab= "Combined activity (log2)",
       xaxt= "n")
  # Legend
  axis(1, padj= -1.25)
  abline(0, 1, lty= "11")
  lines(mult$add, mult$mult)
  leg(cdition, .N)
  
  # Residuals
  par(mai= c(.9,1.25,.9,1.35),
      pty= "m")
  vl_boxplot(residuals_wt,
             residuals_mut,
             compute.pval= list(c(1,2)),
             ylab= "Residuals (log2)",
             xaxt= "n",
             notch= T,
             col= Cc)
  
  # Individual and combined sctivities
  par(mai= c(.9,.9,.9,.9))
  box <- vl_boxplot(unique(.SD[, .(L, indL_wt)])[[2]],
                    unique(.SD[, .(L, indL_mut)])[[2]],
                    unique(.SD[, .(R, indR_wt)])[[2]],
                    unique(.SD[, .(R, indR_mut)])[[2]],
                    log2FoldChange_wt,
                    log2FoldChange_mut,
                    compute.pval= list(c(1,2), c(3,4), c(5,6)),
                    ylab= "Activity (log2)",
                    xaxt= "n",
                    notch= T,
                    at= c(1,2,4,5,7,8),
                    col= Cc)
  axis(1, 
       c(1.5, 4.5, 7.5),
       labels = paste0(c("5'", "3'", "Pairs"), "\nn=", lengths(box$dat1)), 
       line= -.15,
       lwd= 0)
  .SD
}, cdition]
dev.off()