setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(data.table)
source("git_peSTARRSeq/functions/plot_transgene.R")

#---------------------------------------------#
# Import full data
#---------------------------------------------#
meta <- fread("Rdata/metadata_processed.txt")
unique(meta[, .(vllib, spacer, CP, library)])
meta <- meta[vllib %in% c("vllib015", "vllib018", "vllib023",
                          "vllib016", "vllib020", "vllib024") & (DESeq2)]
dat <- unique(meta[, .(CP, spacer, FC_file_ratio)])
dat[, spacer:= switch(spacer,
                      "SCR1"= "no_intron",
                      "shortened-intron4"= "short_intron",
                      "intron4"= "long_intron"), spacer]
dat <- dat[, fread(FC_file_ratio), (dat)]

#---------------------------------------------#
# Import scaled data (smallest sequencing)
#---------------------------------------------#
# dat <- readRDS("Rdata/final_results_table_spacer_size.rds")
# dat[, spacer:= switch(spacer,
#                       "SCR1"= "no_intron",
#                       "shortened-intron4"= "short_intron",
#                       "intron4"= "long_intron"), spacer]

#---------------------------------------------#
# Format
#---------------------------------------------#
dat <- dat[grepl("^dev|^hk", L)
           & grepl("^dev|^hk", R)]
dat[, diff:= log2FoldChange-additive]
dat[, class:= factor(
  fcase(grepl("dev", L) & grepl("dev", R), "dev/dev",
        grepl("dev", L) & grepl("hk", R), "dev/hk",
        grepl("hk", L) & grepl("dev", R), "hk/dev",
        grepl("hk", L) & grepl("hk", R), "hk/hk"),
  levels= c("hk/hk",
            "hk/dev",
            "dev/hk",
            "dev/dev"))]
dat[class=="dev/dev", Cc:= "#74C27A"]
dat[class=="hk/dev", Cc:= "royalblue2"]
dat[class=="dev/hk", Cc:= "cyan"]
dat[class=="hk/hk", Cc:= "tomato"]
dat[, CP:= factor(CP, c("DSCP", "RpS12"))]
dat[, spacer:= factor(spacer, c("no_intron", "short_intron", "long_intron"))]

#---------------------------------------------#
# Plot
#---------------------------------------------#
pdf("pdf/draft/intron_spacers.pdf",
    height= 3,
    width = 3.75)
layout(matrix(1:8, 
              nrow= 2, 
              byrow = T),
       widths = c(1.3,1,1,1),
       heights = c(1, 1.375))
for(var in c("log2FoldChange", "diff"))
{
  ylab <- NA
  if(var=="log2FoldChange")
  {
    ylim <- c(-5, 15)
    ylab <- "Activity (log2)"
  }else if(var=="diff")
  {
    ylim <- c(-7, 6)
    ylab <- "Obs./Exp. Add. (log2)"
  }
  dat[,
      {
        margins <- c(1,2,1.5,0.1)
        xaxt <- "n"
        if(class=="hk/hk")
          margins[2] <- 4
        if(CP=="RpS12")
        {
          margins[c(1, 3)] <- c(6, 0.05)
          xaxt <- "o"
        }
        par(mar= margins,
            lwd= 0.5,
            mgp= c(1.5, 0.5, 0),
            tcl= -0.2,
            las= 1,
            cex= 0.6)
        vl_boxplot(split(get(var), spacer),
                   ylim= ylim,
                   boxcol= Cc,
                   xaxt= xaxt,
                   compute_pval= list(c(1,2), c(2,3)),
                   tilt.names= T,
                   notch= T,
                   ylab= ifelse(class=="hk/hk", ylab, NA))
        transgene(grconvertX(1, "lines", "user"),
                  -5, 
                  CP= ifelse(CP=="DSCP", "dev", "hk"), 
                  rotate = T,
                  cex= 0.5)
        abline(h= 0,
               lty= 2)
        if(CP=="DSCP")
          title(class)
        print("")
      }, keyby= .(CP, class, Cc)]
}
dev.off()

