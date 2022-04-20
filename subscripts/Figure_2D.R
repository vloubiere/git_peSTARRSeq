setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")

pdf("pdf/draft/Figure_2D.pdf",
    width= 6.5,
    height = 5)
par(mgp= c(1.5,0.35,0))
for(side in c("mot_enr_L", "mot_enr_R"))
{
  if(side=="mot_enr_L")
    par(mar= c(7,17,2,6)) else if(side=="mot_enr_R")
      par(mar= c(3,17,2,6))
  # Extract data
  enr <- dat[[side]]$enr[padj<fcase(side=="mot_enr_L", 0.001,
                                    side=="mot_enr_R", 0.01)]
  enr[vl_Dmel_motifs_DB_full, 
      c("motif_cluster_name", "motif_name", "pwm"):= .(i.Motif_cluster_name, i.motif_name, pwms_perc),
      on= "variable==uniqName_noSpecialChar"]
  enr <- enr[, .SD[which.max(abs(log2OR))], motif_cluster_name]
  setnames(enr, 
           c("variable", "motif_cluster_name"), 
           c("uniq_name", "variable"))
  setorderv(enr, "log2OR", -1)
  enr[, check:= fcase(log2OR>0, .I<=10, 
                      log2OR<0, .N-.I<=10)]
  enr <- enr[(check), !"check"]
  
  # Plot
  bar <- plot(enr, 
              axes= F, 
              xlab= "Strong / Weak odd Ratio (log2)", 
              xlim= c(-2.5,2.5), 
              col= c("#00AEEF", "#EE2A7B"))
  axis(1, 
       at= c(-1,0,1),
       labels = c(-1,0,1), 
       tck=-0.02)
  abline(v= 0)
  enr[bar, bar:= i.bar, on= "variable"]
  enr[, top:= bar+0.4]
  enr[, width:= strwidth("M")*0.8*ncol(as.matrix(pwm[[1]])), variable]
  enr[, right:= par("usr")[1]-max(strwidth(paste0("M", variable)))]
  enr[, left:= right-width, variable]
  enr[, height:= 0.8]
  enr[, {
    vl_seqlogo(as.matrix(pwm[[1]]),
               xleft = left,
               ytop = top,
               width = width,
               height= height)
    segments(left,
             top-height,
             left+width,
             top-height,
             lwd= 0.25,
             xpd= T,
             col= "grey")
  }, .(variable, top, left, width, height)]
  left <- min(enr$left)-strwidth("M")*1.5
  top <- enr[, max(top), sign(log2OR)]$V1
  width <- strwidth("M")*1.3
  height <- enr[, diff(range(bar)), sign(log2OR)]$V1
  rect(c(-2, 0),
       par("usr")[4],
       c(0, 2),
       par("usr")[4]+strheight("M")*2.5,
       xpd= T,
       border= NA,
       col= adjustcolor(c("grey70", "grey20"), 0.7))
  text(c(-1, 1),
       par("usr")[4]+strheight("M")*1.25,
       c("Enriched\nin Weak", 
         "Enriched\nin Strong"),
       xpd= T,
       cex= 0.8)
}
dev.off()