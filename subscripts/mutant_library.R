FC <- readRDS("db/FC_tables/vllib029_DESeq2.rds")
FC[, c("group.L", "mut.L", "ID.L"):= tstrsplit(L, "_", keep= c(1,2,4))]
FC[, c("group.R", "mut.R", "ID.R"):= tstrsplit(R, "_", keep= c(1,2,4))]
FC[, ID.L:= paste0(group.L, "_", ID.L)]
FC[, ID.R:= paste0(group.R, "_", ID.R)]
.lm <- lm(log2FoldChange~indL*indR, 
          FC[mut.L=="WT" &  mut.R=="WT"])
FC[, pred:= predict(.lm, .SD)]
FC[, res:= log2FoldChange-pred]


WT <- FC[mut.L=="WT" & mut.R=="WT"]
mut <- FC[!(mut.L=="WT" & mut.R=="WT")]
dat <- merge(WT, 
             mut,
             by= c("ID.L", "ID.R"),
             suffixes= c(".WT", ".mut"))
cols <- c("mut.L.mut", "mut.R.mut")
dat[, (cols):= lapply(.SD, function(x) gsub("add2|add3", "add", x)), .SDcols= cols]
dat <- dat[mut.L.mut=="WT" | mut.R.mut=="WT" 
           | (mut.L.mut %in% c("addTrl", "addTwist") & mut.R.mut %in% c("addTrl", "addTwist"))
           | (mut.L.mut %in% c("mutTrl", "mutTwist") & mut.R.mut %in% c("mutTrl", "mutTwist"))
           | (mut.L.mut=="addDref" & mut.R.mut=="addDref")
           | (mut.L.mut=="mutDref" & mut.R.mut=="mutDref")]


libRes <- FC[mut.L=="WT" & mut.R=="WT" & group.L==group.R & !(group.L %in% c("DHS", "control"))]
libRes[, group:= switch(group.L, 
                        "noMotifAct"= "Trl- Twist- enh.",
                        "shared"= "Shared enh.",
                        "Trl"= "Trl enh.",
                        "twist"= "Twist enh."), group.L]
libRes[, group:= factor(group, 
                        c("Trl- Twist- enh.",
                          "Trl enh.",
                          "Twist enh.",
                          "Shared enh."))]
pdf("pdf/draft/residual_mutant_library.pdf", 
    height = 2.5, 
    width = 5)
par(mar= c(5,4,2,1),
    las= 1,
    tcl= -0.2,
    mgp= c(2, 0.5, 0),
    mfrow= c(1,4))
vl_boxplot(indL~group,
           libRes,
           tilt.names= T, 
           notch= T, 
           compute_pval= list(c(1,2), c(1,3), c(1,4)),
           ylab= "Individual 5' (log2)",
           col= "lightgrey")
vl_boxplot(indR~group,
           libRes,
           tilt.names= T, 
           notch= T, 
           compute_pval= list(c(1,2), c(1,3), c(1,4)),
           ylab= "Individual 3' (log2)",
           col= "lightgrey")
vl_boxplot(log2FoldChange~group,
           libRes,
           tilt.names= T, 
           notch= T, 
           compute_pval= list(c(1,2), c(1,3), c(1,4)),
           ylab= "Activity (log2)",
           col= "lightgrey")
abline(h= median(libRes[group=="Trl- Twist- enh.", log2FoldChange]))
vl_boxplot(res~group,
           libRes,
           tilt.names= T, 
           notch= T, 
           compute_pval= list(c(1,2), c(1,3), c(1,4)),
           ylab= "Observed-Expected (log2)",
           col= "lightgrey")
abline(h= median(libRes[group=="Trl- Twist- enh.", res]))
dev.off()

pdf("pdf/draft/mutant_library.pdf", 4.5, 4.5)
par(mar= c(8,3,4,9),
    las= 1,
    tcl= -0.2,
    mgp= c(1.5, 0.5, 0))
# dat[actL.WT!="Inactive" & actR.WT!="Inactive", {
dat[, {
  .c <- na.omit(.SD[, .(ID.L, ID.R, indL.WT, indL.mut, indR.WT, indR.mut, log2FoldChange.WT, log2FoldChange.mut)])
  if(nrow(.c)>10)
  {
    sub <- .c[, indL.WT:log2FoldChange.mut]
    vl_boxplot(sub, 
               compute_pval = list(c(1,2), c(3,4), c(5,6)),
               tilt.names = T,
               violin= T,
               ylab= "Activity (log2)",
               xaxt= "n", 
               at= c(1.25, 1.75, 3.25, 3.75, 5.25, 5.75), 
               viocol= c("lightgrey", "tomato"),
               # ylim= c(0, max(unlist(sub))),
               col= "white")
    axis(1, 
         c(1.5, 3.5, 5.5), 
         c(paste0("5'\n(", length(unique(ID.L)), ")"), 
           paste0("3'\n(", length(unique(ID.R)), ")"), 
           paste0("Pairs\n(", nrow(.c), ")")),
         lty= 0)
    legend(par("usr")[2],
           par("usr")[4],
           fill= c("lightgrey", "tomato"),
           legend= c("WT x WT", 
                     paste0(mut.L.mut, " x ", mut.R.mut)),
           bty= "n",
           xpd= T,
           cex= 0.8)
  }
  .SD
}, .(mut.L.mut, mut.R.mut)]
dev.off()

file.show("pdf/draft/mutant_library.pdf")

pdf("pdf/draft/mutant_library_residuals.pdf", 4, 4.5)
par(mar= c(8,8,4,8),
    las= 1,
    tcl= -0.2,
    mgp= c(2.5, 0.5, 0))
dat[, {
  vl_boxplot(res.WT, 
             res.mut,
             compute_pval= list(c(1,2)),
             notch= T,
             xaxt= "n",
             col= c("lightgrey", "tomato"),
             ylab= "Observed-Expected (log2)")
  text(1.5, 
       par("usr")[3]-strheight("M"), 
       paste0("n= ", sum(!is.na(res.mut))),
       xpd= T)
  legend(par("usr")[2],
         par("usr")[4],
         fill= c("lightgrey", "tomato"),
         legend= c("WT x WT", 
                   paste0(mut.L.mut, " x ", mut.R.mut)),
         bty= "n",
         xpd= T,
         cex= 0.8)
  .SD
}, .(mut.L.mut, mut.R.mut)]
dev.off()

file.show("pdf/draft/mutant_library_residuals.pdf")