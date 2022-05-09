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
dat <- lib[vllib %in% c("vllib015", "vllib028")][, diff:= log2FoldChange-additive]
dat <- merge(dat[vllib=="vllib015"],
             dat[vllib=="vllib028"],
             by= c("L", "R"),
             suffixes= c("_low", "_high"))
dat <- dat[class_low=="enh./enh."]

pdf("pdf/draft/Figure_4B.pdf",
    height= 3,
    width = 5)
par(mfrow= c(1,3),
    mar= c(3,5,1,0),
    mgp= c(1.5, 0.5, 0),
    las= 1,
    tcl= -0.2)
vl_boxplot(unique(dat[, .(L, 
                          low= median_L_low, 
                          high= median_L_high)][, !"L"]),
           compute_pval= list(c(1,2)),
           tilt.names= T,
           ylab= "Individual activity (log2)",
           main= "5' enhancers")
text(x= par("usr")[1],
     y= mean(c(par("usr")[3], grconvertY(0, "nfc", "user"))),
     "CP basal\nactivity",
     pos= 2,
     xpd= T)
vl_boxplot(unique(dat[, .(R, 
                          low= median_R_low, 
                          high= median_R_high)][, !"R"]),
           compute_pval= list(c(1,2)),
           tilt.names= T,
           ylab= "Individual activity (log2)",
           main= "3' enhancers")
vl_boxplot(dat[, .(low= diff_low, 
                   high= diff_high)], 
           compute_pval= list(c(1,2)), 
           tilt.names= T,
           ylab= "obs./Exp. Add. (log2)",
           main= "Enhancer pairs")
dev.off()
