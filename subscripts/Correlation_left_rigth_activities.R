setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import data ----
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
dat <- merge(unique(dat[, .(ID= L, indL)]),
             unique(dat[, .(ID= R, indR)]))

# Add TWIST data ----
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
twist <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/BA_300bp_TWIST_STARRSeq.txt")
lib[twist, dev_log2FC_TWIST:= i.dev_log2FoldChange, on= "BA_ID==ID"]
dat[lib, dev_log2FC_TWIST:= i.dev_log2FC_TWIST, on= "ID==ID_vl"]
dat[, TWIST_class:= cut(dev_log2FC_TWIST, 
                        c(-Inf, 2, 6, 8, Inf),
                        c("Inactive", "Weak", "Medium", "Strong"))]
dat[grepl("^control", ID) | is.na(TWIST_class), TWIST_class:= ifelse(grepl("^control", ID), "Control", "Inactive")]
dat[, TWIST_class:= factor(TWIST_class, 
                           c("Control", "Inactive", "Weak", "Medium", "Strong"))]
dat[, TWIST_col:= c("grey20", "#4D9221", "#B8E186", "#F1B6DA", "#C51B7D")[TWIST_class]]
setorderv(dat, "TWIST_class", order= -1)

# Model ----
.lm <-lm(indL~indR, dat)

# Plot ----
pdf("pdf/draft/Correlation_left_rigth_activities.pdf", 
    width = 3, 
    height = 2.9)
# Scatter plot
par(mar= c(3.5, 3, 0.5, 0.5), 
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2,
    bty= "n")
dat[, {
  # Scatter plot
  plot(indR,
       indL,
       pch= 16,
       cex= 0.5,
       col= adjustcolor(TWIST_col, 0.7),
       las= 1,
       xlab= "3' activity (log2)",
       ylab= "5' activity (log2)")
}]
# Legend
abline(.lm, lty= "11")
unique(dat[, .(TWIST_class, TWIST_col)])[,{
  legend("bottomright",
         col= TWIST_col,
         pch= 16,
         legend= TWIST_class,
         bty= "n",
         cex= 0.7)
}]
vl_plot_R2(rsquare = summary(.lm)$r.squared,
           inset= c(-0.1,0))
dev.off()