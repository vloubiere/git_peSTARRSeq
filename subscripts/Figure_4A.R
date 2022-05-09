setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(data.table)
source("git_peSTARRSeq/functions/plot_transgene.R")

#---------------------------------------------#
# Import full data
#---------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")

#---------------------------------------------#
# Format
#---------------------------------------------#
dat <- lib[vllib %in% c("vllib002", "vllib006")][, diff:= log2FoldChange-additive]
dat <- merge(dat[vllib=="vllib002"],
             dat[vllib=="vllib006"],
             by= c("L", "R"),
             suffixes= c("_short", "_long"))
dat <- dat[class_short=="enh./enh." & class_long=="enh./enh."]

pdf("pdf/draft/Figure_4A.pdf",
    height= 3,
    width = 5)
par(mfrow= c(1,3),
    mar= c(3,5,1,0),
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2)
vl_boxplot(unique(dat[, .(L, 
                          `300bp`= median_L_short, 
                          `2kb`= median_L_long)][, !"L"]),
           compute_pval= list(c(1,2)),
           tilt.names= T,
           ylab= "Individual activity (log2)",
           main= "5' enhancers")
text(x= par("usr")[1],
     y= mean(c(par("usr")[3], grconvertY(0, "nfc", "user"))),
     "Spacer size",
     pos= 2,
     xpd= T)
vl_boxplot(unique(dat[, .(R, 
                          `300bp`= median_R_short, 
                          `2kb`= median_R_long)][, !"R"]),
           compute_pval= list(c(1,2)),
           tilt.names= T,
           ylab= "Individual activity (log2)",
           main= "3' enhancers")
vl_boxplot(dat[, .(`300bp`= diff_short, 
                   `2kb`= diff_long)], 
           compute_pval= list(c(1,2)), 
           tilt.names= T,
           ylab= "obs./Exp. Add. (log2)",
           main= "Enhancer pairs")
dev.off()