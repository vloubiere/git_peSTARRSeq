setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("db/FC_tables_DESeq2/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_DESeq2_final_oe.rds")
# Retrieve candidates vs control sequences
pl <- rbind(lib[, .(ID= L, ind= indL, actClass= actClassL, side= "5'")],
            lib[, .(ID= R, ind= indR, actClass= actClassR, side= "3'")])
pl <- unique(pl)
pl[, group:= factor(ifelse(grepl("^control", ID), "rdm seq.", "Candidates"),
                           c("rdm seq.", "Candidates"))]
# Split Active enhancers based on their strength
pl[, act:= fcase(actClass=="inactive", "Inactive",
                 ind<=2, "Low",
                 ind<=4, "Medium",
                 default= "Strong")]
pl[, act:= factor(act, c("Inactive", "Low", "Medium", "Strong"))]
pl <- pl[, {tot <- .N; .SD[, .(.N, perc= .N/tot*100), act]}, .(side, group)]
pl[, ypos:= cumsum(perc)-perc/2, .(side, group)]

Cc <- adjustcolor(c("cornflowerblue", "gold", "tomato", "red"), 0.5)


pdf("pdf/draft/Candidates_classification_vllib002.pdf",
    height = 3.25, 
    width = 4)
# Scatter plot
par(mar= c(4, 3, 1.5, 1),
    mfrow= c(1,2),
    oma= c(0,0,0,3),
    mgp= c(2, 0.5, 0),
    las= 2,
    tcl= -0.2,
    lwd= 0.25)
pl[, {
  bar <- barplot(perc~act+group,
                 .SD,
                 col= Cc,
                 xlab= NA,
                 main= side,
                 ylab= "% of oligos",
                 xaxt= "n")
  vl_tilt_xaxis(bar, labels= levels(group))
  text(bar[group],
       ypos,
       ifelse(N>10, N, NA),
       cex= 0.7)
  if(.GRP==.NGRP)
    legend(par("usr")[2],
           par("usr")[4],
           legend = rev(levels(act)),
           fill= rev(Cc),
           xpd= NA,
           bty= "n",
           cex= 0.7,
           border= NA)
  ""
}, side]
dev.off()
