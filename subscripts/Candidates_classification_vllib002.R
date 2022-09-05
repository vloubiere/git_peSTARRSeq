setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]

# Retrieve candidates vs control sequences
pl <- data.table(enh= unique(lib[, c(L, R)]))
pl[, group:= fcase(grepl("^control", enh), "rdm seq.",
                   default = "Candidates")]
pl[, group:= factor(group, c("rdm seq.", "Candidates"))]

# Merge individual act
pl <- rbindlist(list("5'"= pl[unique(lib[, .(L, class= class_act_L, act= median_L)]), on= "enh==L"],
                     "3'"= pl[unique(lib[, .(R, class= class_act_R, act= median_R)]), on= "enh==R"]), 
                idcol = T)
pl[, .id:= factor(.id, c("5'", "3'"))]
# Split Active enhancers based on their strength
pl[, act:= fcase(class=="inactive", "inactive",
                 act<=2, "low (<=2)",
                 act<=4, "medium (2-4)",
                 default= "strong (>4)")]
pl[, act:= factor(act, c("inactive", "low (<=2)", "medium (2-4)", "strong (>4)"))]

pdf("git_peSTARRSeq/subscripts/Candidates_classification_vllib002.pdf",
    width = 3, 
    height = 3)
# Scatter plot
par(mar= c(2.5, 3, 1.5, 5), 
    mgp= c(2, 0.5, 0),
    las= 1,
    tcl= -0.2)
pl[, {
  # Compute number and perc
  tab <- .SD[, {
    tot_group <- .N
    .SD[, .(count= .N, perc= .N/tot_group*100), act]
  }, .(.id, group)]
  # Barplot percentage
  perc <- dcast(tab, 
               act~.id+group, 
               value.var = "perc")
  perc <- as.matrix(perc, 1)
  Cc <- adjustcolor(c("cornflowerblue", "gold", "tomato", "red"), 0.5)
  bar <- barplot(perc, 
                 col= Cc,
                 names.arg = rep(NA, ncol(perc)),
                 space= c(1,0.1),
                 ylab= "% of oligos", 
                 border=NA)
  text(bar,
       par("usr")[3]-strheight("M", cex= 0.25), 
       labels = gsub("5'_|3'_", "", colnames(perc)),
       pos= 2,
       srt= 45,
       xpd= T,
       offset= -0.25)
  # Add numbers
  count <- dcast(tab, 
                act~.id+group, 
                value.var = "count")
  count <- as.matrix(count, 1)
  cumPerc <- apply(perc, 2, cumsum)-perc/2
  text(rep(bar, each= nrow(perc)),
       c(cumPerc), 
       labels = ifelse(c(perc>5), c(count), NA),
       xpd= T,
       cex= 0.7)
  # Legend
  segments(bar[c(1,3)],
           par("usr")[4]+strheight("M", cex= 0.5),
           bar[c(2,4)],
           par("usr")[4]+strheight("M", cex= 0.5),
           xpd= T)
  text(rowMeans(matrix(bar, byrow = T, ncol= 2)),
       par("usr")[4]+strheight("M"),
       labels= levels(tab$.id),
       xpd= T,
       pos= 3,
       offset= 0)
  legend(par("usr")[2],
         par("usr")[4],
         fill= rev(Cc), 
         legend = rev(levels(tab$act)),
         xpd= T,
         bty= "n",
         cex= 0.7,
         border= NA)
}]
dev.off()
