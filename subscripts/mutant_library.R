setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_mutant_library_FC_DESeq2.rds")

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("groupL", "mutL", "L"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "R"):= tstrsplit(R, "_", keep= c(1,2,4))]
dat[, L:= paste0(groupL, L)]
dat[, R:= paste0(groupR, R)]

# Compute expected multiplicative ----
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
dat[, `linear model`:= predict(model, .SD)]
dat[, residuals:= log2FoldChange-`linear model`]

# Merge WT and mutant variants ----
pl <- merge(dat[mutL=="WT" & mutR=="WT", .(L, R, indL, indR, log2FoldChange, `linear model`, residuals)],
            dat[mutL!="WT" & mutR!="WT", .(L, R, mutL, mutR, indL, indR, log2FoldChange, `linear model`, residuals)],
            by= c("L", "R"),
            suffixes= c("_wt", "_mut"))

# Retrieve combinations of interest ----
# Pasting Trl/Twist motifs in enhancers with no motifs
pl[grepl("^noMotifAct", L) & grepl("^noMotifAct", R), 
   cdition:= fcase(grepl("add.*Trl", mutL) & grepl("add.*Trl", mutR), "Added Trl motifs",
                   grepl("add.*Twist", mutL) & grepl("add.*Twist", mutR), "Added Twist motifs")]
# Pasting Dref motif in dev enhancers
pl[grepl("^noMotifAct", L) & grepl("^noMotifAct", R) & grepl("add.*Dref", mutL) & grepl("add.*Dref", mutR), cdition:= "Added Dref motifs"]
# Mutate motifs in enhancers that contain them
pl[grepl("mut.*Trl", mutL) & grepl("mut.*Trl", mutR), cdition:= "Mutated Trl motifs",]
pl[grepl("mut.*Twist", mutL) & grepl("mut.*Twist", mutR), cdition:= "Mutated Twist motifs"]
# pl[grepl("mut.*Dref", mutL) & grepl("mut.*Dref", mutR), cdition:= "Mutated Dref motifs"] # NOT USED
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
  xlim <- range(c(`linear model_wt`, `linear model_mut`))
  ylim <- range(c(log2FoldChange_wt, log2FoldChange_mut))
  wt <- unique(.SD[, .(indL_wt,
                       indR_wt,
                       `linear model_wt`,
                       log2FoldChange_wt)])
  wt[, {
    plot(`linear model_wt`,
         log2FoldChange_wt, 
         main= cdition,
         cex= .5,
         pch= 16,
         col= Cc[1],
         xlim= xlim,
         ylim= ylim,
         xlab= "Predicted (log2)",
         ylab= "Activity (log2)",
         xaxt= "n")
    axis(1, padj= -1.25)
  }]
  abline(0, 1, lty= "11")
  points(`linear model_mut`,
         log2FoldChange_mut,
         pch= 16,
         col= adjustcolor(Cc[2], .6),
         cex= .5)
  # Legend
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