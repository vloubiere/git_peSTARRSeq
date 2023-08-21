setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/vllib029_DESeq2.rds")
dat[, c("groupL", "mutL", "IDL"):= tstrsplit(L, "_", keep= c(1,2,4))]
dat[, c("groupR", "mutR", "IDR"):= tstrsplit(R, "_", keep= c(1,2,4))]
# Make sequence ID unique (numbers are only unique within one group)
dat[, IDL:= paste0(groupL, IDL)]
dat[, IDR:= paste0(groupR, IDR)]
.lm <- lm(log2FoldChange~indL*indR,
          dat[mutL=="WT" & mutR=="WT"])
dat[, predicted:= predict(.lm, .SD)]
dat[, residuals:= log2FoldChange-predicted]
smoothScatter(dat[mutL=="WT" & mutR=="WT", .(predicted, log2FoldChange)])
abline(0,1)

# Select active pairs
active <- dat[mutL=="WT" & indL>0.5 & mutR=="WT" & indR>0.5, .(IDL, IDR)]
active <- merge(active, dat, by= c("IDL", "IDR"))

# Extract data for combinations of interest ----
cmb <- data.table(pattern= c("add.*Trl", "add.*Twist", "add.*Dref",
                             "mut.*Trl", "mut.*Twist", "mut.*Dref"),
                  cdition= c("Added Trl motifs", "Added Twist motifs", "Added Dref motifs",
                             "Mutated Trl motifs", "Mutated Twist motifs", "Mutated Dref motifs"))
cmb <- cmb[, {
  .c <- dat[grepl(pattern, mutL) & grepl(pattern, mutR)]
  .ctl <- dat[mutL=="WT" & mutR=="WT"]
  res <- merge(.ctl[, .(IDL, IDR, indL, indR, log2FoldChange, predicted, residuals)],
               .c[, .(IDL, IDR, log2FoldChange, predicted, residuals)],
               by= c("IDL", "IDR"),
               suffixes= c("_ctl", "_mut"))
  res <- melt(res,
              measure.vars = patterns("log2FoldChange"= "^log2FoldChange",
                                      "predicted"= "^predicted",
                                      "residuals"= "^residuals"))
  res[, variable:= c("WT", "Mutated")[variable]]
  res
}, cdition]
cmb[, variable:= factor(variable, c("WT", "Mutated"))]


# Plot results ----
cmb <- rbindlist(list(allPairs= cmb,
                      actPairs= cmb[indL>=0.5 & indR>=0.5]),
                 idcol = "group")
cmb[, pdf:= paste0("pdf/draft/mutant_library_act_res_", group, ".pdf")]
Cc <- c("lightgrey", "pink1")
cmb[, {
  pdf(pdf, 4.25, 2.75)
  par(mar= c(2,3,0,0),
      oma= c(0,0,3,9),
      mgp= c(1.5, 0.5, 0),
      mfrow= c(1,2),
      tcl= -0.2,
      las= 1)
  .SD[, {
    vl_boxplot(log2FoldChange~variable,
               compute_pval= list(c(1,2)),
               ylab= "Activity (log2)",
               notch= T,
               xaxt= "n",
               col= Cc)
    vl_boxplot(residuals~variable,
               compute_pval= list(c(1,2)),
               ylab= "Residuals (log2)",
               xaxt= "n",
               notch= T,
               col= Cc)
    text(x= grconvertX(0.5, "ndc", "user"),
         y= grconvertY(1, "ndc", "user"),
         labels= paste0(unique(cdition), " (n=", formatC(.N/2, big.mark = ","), ")"),
         pos= 1,
         xpd= NA)
    legend(par("usr")[2],
           par("usr")[4],
           unique(variable),
           fill= Cc,
           bty= "n", 
           xpd= NA)
    print(pdf)
  }, cdition]
  dev.off()
}, pdf]