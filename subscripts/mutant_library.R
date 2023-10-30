setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/vllib029_DESeq2.rds")

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("groupL", "mutL", "IDL"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "IDR"):= tstrsplit(R, "_", keep= c(1,2,4))]
dat[, IDL:= paste0(groupL, IDL)]
dat[, IDR:= paste0(groupR, IDR)]

# Compute expected multiplicative ----
model <- readRDS("db/linear_models/lm_vllib002.rds")
dat[, multiplicative:= predict(model, .SD)]
dat[, residuals:= log2FoldChange-multiplicative]

# Extract data for combinations of interest ----
pl <- merge(dat[mutL=="WT" & mutR=="WT", .(IDL, IDR, indL, indR, log2FoldChange, residuals)],
            dat[mutL!="WT" & mutR!="WT", .(IDL, IDR, mutL, mutR, indL, indR, log2FoldChange, residuals)],
            by= c("IDL", "IDR"),
            suffixes= c("_wt", "_mut"))
pl[, cdition:= fcase(grepl("add.*Trl", mutL) & grepl("add.*Trl", mutR), "Added Trl motifs",
                     grepl("add.*Twist", mutL) & grepl("add.*Twist", mutR), "Added Twist motifs",
                     grepl("add.*Dref", mutL) & grepl("add.*Dref", mutR), "Added Dref motifs",
                     grepl("mut.*Trl", mutL) & grepl("mut.*Trl", mutR), "Mutated Trl motifs",
                     grepl("mut.*Twist", mutL) & grepl("mut.*Twist", mutR), "Mutated Twist motifs",
                     grepl("mut.*Dref", mutL) & grepl("mut.*Dref", mutR), "Mutated Dref motifs",
                     default= "mixed")]

# Only consider enhancers that were active in the first place ----
pl <- pl[(grepl("^noMotifAct", IDL) & grepl("^noMotifAct", IDL)) 
         | (grepl("^Trl", IDL) & grepl("^Trl", IDL))
         | (grepl("^twist", IDL) & grepl("^twist", IDL))
         | (grepl("shared", IDL) & grepl("shared", IDR))]

pdf("pdf/draft/mutant_library_per_cdition.pdf", 
    width = 6,
    height= 3)
par(mfrow= c(1,2),
    mgp= c(1, 0.25, 0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    font.main= 1,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75)
Cc <- c("lightgrey", "pink1")
pl[cdition!="mixed", {
  # Residuals
  par(mai= c(.9,1.25,.9,1.35))
  vl_boxplot(residuals_wt,
             residuals_mut,
             compute_pval= list(c(1,2)),
             ylab= "Residuals (log2)",
             xaxt= "n",
             notch= T,
             col= Cc)
  leg <- substitute(legend(par("usr")[2],
                           par("usr")[4],
                           fill= Cc,
                           legend= c("WT",
                                     fcase(grepl("Mutated", cdition), "Mutant",
                                           grepl("Added", cdition), "Added motif")),,
                           bty= "n",
                           xpd= NA,
                           cex= 7/12))
  eval(leg)
  # Activities
  par(mai= c(.9,.9,.9,.9))
  box <- vl_boxplot(unique(.SD[, .(IDL, indL_wt)])[[2]],
                    unique(.SD[, .(IDL, indL_mut)])[[2]],
                    unique(.SD[, .(IDR, indR_wt)])[[2]],
                    unique(.SD[, .(IDR, indR_mut)])[[2]],
                    log2FoldChange_wt,
                    log2FoldChange_mut,
                    compute_pval= list(c(1,2), c(3,4), c(5,6)),
                    ylab= "Activity (log2)",
                    xaxt= "n",
                    notch= T,
                    at= c(1,2,4,5,7,8),
                    col= Cc)
  text(x= grconvertX(0.5, "ndc", "user"),
       y= grconvertY(1, "ndc", "user"),
       labels= paste0(unique(cdition)," (all, n=", formatC(.N, big.mark = ","), ")"),
       pos= 1,
       xpd= NA)
  axis(1, 
       c(1.5, 4.5, 7.5),
       labels = paste0(c("5'", "3'", "Pairs"), "\nn=", lengths(box$dat1)), 
       line= -.15,
       lwd= 0)
  eval(leg)
  .SD
}, cdition]
dev.off()