setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
# Retrieve candidates vs control sequences
pl <- rbind(lib[, .(ID= L, group= groupL, act= actL, side= "5'")],
            lib[, .(ID= R, group= groupR, act= actR, side= "3'")])
pl <- unique(pl)
pl[, group:= factor(group, c("Random seq.", "Candidate seq."))]
# Split Active enhancers based on their strength
pl[, act:= factor(act, c("Inactive", "Low", "Medium", "Strong"))]
pl <- pl[, {tot <- .N; .SD[, .(.N, perc= .N/tot*100), act]}, .(side, group)]
pl[, ypos:= cumsum(perc)-perc/2, .(side, group)]

Cc <- adjustcolor(c("cornflowerblue", "gold", "tomato", "red"), 0.5)

pdf("pdf/draft/Candidates_classification_vllib002.pdf",
    height = 3.5, 
    width = 3.75)
# Scatter plot
par(mar= c(4.75, 1, 1.5, 1),
    mfrow= c(1,2),
    oma= c(0,3,0,3),
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
                 xaxt= "n")
  vl_tilt_xaxis(bar, labels= levels(group))
  text(bar[group],
       ypos,
       ifelse(N>10, N, NA),
       cex= 0.7)
  if(.GRP==1)
    title(ylab= "% of oligos", xpd= NA)else
      legend(par("usr")[2],
             par("usr")[4],
             legend = rev(levels(act)),
             fill= rev(Cc),
             xpd= NA,
             bty= "n",
             cex= 0.7,
             border= NA)
  .SD
}, side]
dev.off()

file.show("pdf/draft/Candidates_classification_vllib002.pdf")
