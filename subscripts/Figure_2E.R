setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")

pdf("pdf/draft/Figure_2E.pdf",
    width= 4,
    height = 3)
par(mgp= c(1.5,0.35,0),
    cex= 0.66)
for(side in c("mot_enr_L", "mot_enr_R"))
{
  if(side=="mot_enr_L")
    par(mar= c(7,14,2,6)) else if(side=="mot_enr_R")
      par(mar= c(2.75,16,1.75,6))
  # Extract data
  enr <- dat[[side]]$enr[padj<fcase(side=="mot_enr_L", 0.001,
                                    side=="mot_enr_R", 0.01)]
  # Extract BA motif cluster names
  enr[vl_Dmel_motifs_DB_full, 
      c("motif_cluster_name", "motif_name", "pwm"):= .(i.Motif_cluster_name, i.motif_name, pwms_perc),
      on= "variable==uniqName_noSpecialChar"]
  # For each cluster, keep top enriched motif
  enr <- enr[, .SD[which.max(abs(log2OR))], motif_cluster_name]
  setnames(enr, 
           c("variable", "motif_cluster_name"), 
           c("uniq_name", "variable"))
  setorderv(enr, "log2OR", -1)
  # Keep top N enrichments
  enr[, check:= fcase(log2OR>0, .I<=10,
                      log2OR<0, (.N-.I+1)<=10)]
  enr <- enr[(check), !"check"]

  ##### Plot #####
  bar <- plot(enr, 
              axes= F, 
              xlab= "A / B odd Ratio (log2)", 
              xlim= c(-3,3), 
              col= c("#00AEEF", "#EE2A7B"))
  axis(1, 
       at= c(-3,0,3),
       labels = c(-3,0,3), 
       tck=-0.02)
  # Plot PWMs
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
  # Plot clusters rectangles
  rleft <- c(par("usr")[1], 0)
  rright <- c(0, par("usr")[2])
  rect(rleft,
       par("usr")[4],
       rright,
       par("usr")[4]+strheight("M")*2.5,
       xpd= T,
       lwd= 0.25,
       col= adjustcolor(c("grey60", "grey90"), 0.7))
  text(rleft[1], 
       par("usr")[4]+strheight("M")*1.25,
       switch(side, 
              "mot_enr_L"= "5' enhancer",
              "mot_enr_R"= "3' enhancer"),
       pos= 2,
       xpd= T,
       cex= 2)
  abline(v= 0, 
         lwd= 0.25)
  text(sapply(list(rleft, rright), mean),
       par("usr")[4]+strheight("M")*1.25,
       c("Cluster B", 
         "Cluster A"),
       xpd= T)
}
dev.off()