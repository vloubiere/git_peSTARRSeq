setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")
dat <- dat[!grepl("control", L) & !grepl("control", R) 
           & actL!="Inactive"
           & actR!="Inactive"]

# Select variables of interest
top <- data.table(name= c("Ebox/CATATG/twi", "Trl/1", "DRE/1"),
                  variable= c("flyfactorsurvey__CG16778_SANGER_5_FBgn0003715", 
                              "cisbp__M2014", 
                              "homer__AVYTATCGATAD_DREF"),
                  min_count= c(2,3,1))
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
lib <- lib[ID_vl %in% dat[, c(L, R)]]

# Plot object
pl <- top[, {
  .c <- vl_motif_counts(lib$enh_sequence, 
                        variable, 
                        genome = "dm3", 
                        collapse_overlapping = T)
  enr <- lib$ID_vl[.c[[1]]>=min_count]
  enr <- dat[L %in% enr & R %in% enr]
  print(nrow(enr))
  ctls <- lib$ID_vl[.c[[1]]==0]
  ctls <- dat[L %in% ctls & R %in% ctls]
  ctls <- enr[, {
    idx <- seq(.N)
    res <- .SD[, {
      ctls[, dist:= sqrt((cL-indL)^2+(cR-indR)^2)]
      ctls <- ctls[order(dist)]
      .c <- ctls[1]
      ctls <- ctls[-1]
      .c
    }, .(cL= indL, cR= indR, idx)]
  }]
  print(".")
  cbind(enr[,.(L, R, indL, indR, log2FoldChange)],
        ctls[,.(ctl_indL= indL, ctl_indR= indR, ctl_log2FoldChange= log2FoldChange)])
}, .(variable, name, min_count)]
pl[, name:= factor(name, c("DRE/1", "Trl/1", "Ebox/CATATG/twi"))]

# PLOT
pdf("pdf/draft/activity_predicitve_motifs_activity_matched.pdf", 
    height= 3.5,
    width= 7.5)
par(mar= c(2.1, 7, 2.1, 18),
    mgp= c(1.5,0.5,0),
    oma= c(0,0,0,5),
    las= 1,
    tcl= -0.2,
    xpd= NA)
pl[, {
  xpos <- seq(1, 5, length.out= 3)
  vl_boxplot(ctl_indL,
             indL,
             ctl_indR,
             indR,
             ctl_log2FoldChange,
             log2FoldChange,
             at= rep(xpos, each= 2)+c(-0.35, 0.35),
             tilt.names= T, 
             compute_pval= list(c(1,2), c(3,4), c(5,6)),
             main= paste0(name, " (", 
                          length(unique(L)), " x 5'; ",
                          length(unique(R)), " x 3'; ",
                          .N, " pairs)"),
             xaxt= "n",
             col= c("lightgrey", "rosybrown1"),
             ylab= "Activity (log2)",
             notch= T)
  axis(1, at= xpos, labels = c("5'", "3'", "Pair"))
  legend(par("usr")[2],
         par("usr")[4],
         fill= c("lightgrey", "rosybrown1"),
         legend= c("Controls (no motifs)", 
                   paste0("Motif pairs (>= ", min_count, ")")),
         bty= "n")
  .SD
}, .(name, min_count)]

# Unique plot
par(mar= c(4.1, 18, 1.5, 10),
    mgp= c(1.5,0.5,0),
    oma= c(0,0,0,0))
.c <- unlist(split(pl[, .(ctl_log2FoldChange, log2FoldChange)], pl$name), recursive = F)
vl_boxplot(.c,
           horizontal= T,
           at= rep(seq(levels(pl$name)), each= 2)+c(-0.15, 0.15),
           boxwex = .2,
           col= c("lightgrey", "rosybrown1"),
           compute_pval= list(c(1,2), c(3,4), c(5,6)),
           yaxt= "n",
           xlab= "Pairs' activity (log2)",
           notch= T)
names.arg <- table(pl$name)
names.arg <- paste0(names(names.arg), " (", names.arg, ")")
axis(2, seq(levels(pl$name)), names.arg)
pwms <- vl_Dmel_motifs_DB_full[match((unique(pl[, variable[order(name)]])), motif_ID)]$pwms_perc
vl_seqlogo(lapply(pwms, as.matrix),
           x= par("usr")[1]-max(strwidth(paste0(names.arg, "M"))),
           y= seq(levels(pl$name)))
legend(par("usr")[2],
       par("usr")[4],
       fill= c("lightgrey", "rosybrown1"),
       legend= c("Controls (no motifs)",
                 "Motif pairs"),
       bty= "n")
dev.off()

file.show("pdf/draft/activity_predicitve_motifs_activity_matched.pdf")
