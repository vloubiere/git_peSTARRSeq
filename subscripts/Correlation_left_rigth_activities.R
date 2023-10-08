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
.lm <- lm(indL~indR, dat)

# Plot ----
pdf("pdf/draft/Correlation_left_rigth_activities.pdf", 
    width = 3, 
    height = 3)
par(mai= c(0.75,0.75,0.75,0.75), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2)
dat[, {
  # Scatter plot
  plot(indR,
       indL,
       pch= 16,
       cex= 0.5,
       col= adjustcolor(TWIST_col, 0.7),
       las= 1,
       xaxt= "n",
       xlab= "3' activity (log2)",
       ylab= "5' activity (log2)")
  axis(1, padj = -1.25)
  vl_plot_coeff(value = cor.test(indR, indL)$estimate,
                type= "pcc",
                cex= 7/12)
}]
abline(.lm, lty= "11")
# Legend
unique(dat[, .(TWIST_class, TWIST_col)])[,{
  legend("bottomright",
         col= adjustcolor(TWIST_col, 0.7),
         pch= 16,
         legend= TWIST_class,
         bty= "n",
         cex= 6/12,
         y.intersp = .75,
         inset= c(-.3, 0),
         xpd= T)
}]
dev.off()