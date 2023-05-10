setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
dat <- merge(unique(dat[, .(ID= L, indL)]),
             unique(dat[, .(ID= R, indR)]))
# Add TWIST data
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

#------------------------------------------------------#
# PLOT
#------------------------------------------------------#
pdf("pdf/draft/Correlation_left_rigth_activities.pdf", 
    width = 3, 
    height = 3)
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
       xlab= "3' individual activity (log2)",
       ylab= "5' individual activity (log2)")
  # Linear model
  .lm <- lm(indL~indR)
  abline(.lm, lty= "11")
  
  # Legend
  leg <- unique(.SD[, TWIST_col, keyby= TWIST_class])
  leg[, legend("topleft",
               legend = c(paste0("R2= ", round(summary(.lm)$r.squared, 2)),
                          as.character(rev(TWIST_class))),
               pch= c(NA, 
                      rep(16, length(TWIST_class))),
               col= c("black", 
                      rev(TWIST_col)),
               lty= c("11",
                      rep(NA, length(TWIST_class))),
               cex= 0.8,
               seg.len= 0.5,
               bty= "n",
               y.intersp= 0.8)]
}]
dev.off()

file.show("pdf/draft/Correlation_left_rigth_activities.pdf")