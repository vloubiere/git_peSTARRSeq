setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import ----
dat <- rbindlist(list(hkCP= readRDS("db/FC_tables/RpS12_focused_WT_DESeq2.rds"),
                      devCP= readRDS("db/FC_tables/DSCP_focused_WT_DESeq2.rds")),
                 idcol = "CP")
dat[, CP:= factor(CP, c("hkCP", "devCP"))]
dat <- dat[(grepl("^hk|^dev", L) & grepl("^hk|^dev", R))]
dat[, classL:= tstrsplit(L, "_", keep= 1)]
dat[, classR:= tstrsplit(R, "_", keep= 1)]
dat[, class:= paste0(ifelse(classL=="dev", "Dev.", "Hk."), " / ", ifelse(classR=="dev", "Dev.", "Hk."))]
dat[class %in% c("Hk. / Dev.", "Dev. / Hk."), class:= "Hk. / Dev. & Dev. / Hk."]
dat[, class:= factor(class,
                     c("Hk. / Hk.",
                       "Hk. / Dev. & Dev. / Hk.",
                       "Dev. / Dev."))]
dat[, additive:= log2(2^indL+2^indR-1)]
pl <- melt(dat,
           id.vars = c("CP", "class"),
           measure.vars = c("additive", "log2FoldChange"))

pdf("pdf/draft/review_boxplot_hkCP_dev_hk_heterotypic_pairs.pdf",
    width = 3, 
    height = 3)
par(mai= c(.9, .9, .9, .2), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75,
    xpd= NA)
vl_boxplot(value~variable+class+CP,
           pl,
           tilt.names= T,
           col= c("lightgrey", "white"),
           names= function(x) gsub("..hkCP|..devCP|log2FoldChange.|additive.", "", x),
           compute.pval= list(c(1,2), 
                              c(3,4), 
                              c(5,6), 
                              c(7,8), 
                              c(9,10), 
                              c(11,12)))
dev.off()