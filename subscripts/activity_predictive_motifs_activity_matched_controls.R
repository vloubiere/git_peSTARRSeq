setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

if(!file.exists("db/activity_matched_controls/activity_matched_controls.rds"))
{
  # Import data ----
  dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
  dat <- dat[!grepl("^control", L) & !grepl("^control", R)]
  dat$ctlL <- dat$ctlR <- NULL
  setorderv(dat, c("indL", "indR"))
  mot <- readRDS("db/motif_counts/twist008_motif_counts.rds")
  mot <- melt(mot, "ID")
  mot <- data.table(min_counts= 1:2)[, (mot), min_counts]
  # Chose activity matched controls ----
  res <- mot[, {
    # Define enr (motif counts >0) and ctl sets (no motif)
    posIDs <- ID[value>=min_counts]
    enr <- dat[L %in% posIDs & R %in% posIDs]
    negIDs <- ID[value==0]
    ctls <- dat[!L %in% posIDs & !R %in% posIDs]
    # To limit search space, max diff between enr and ctls
    breaks <- 0.1
    enr[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
    ctls[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
    # While some pairs have no control
    enr[, ctlL:= as.character(NA)]
    tot <- nrow(enr)
    while(nrow(ctls) && any(is.na(enr$ctlL)))
    {
      # For each enr pair, define the ctl pair with smallest euclidean distance
      .c <- ctls[enr[is.na(ctlL)], {
        dist <- sqrt((x.indL-i.indL)^2+(x.indR-i.indR)^2)
        idx <- which.min(dist)
        .(L= i.L,
          R= i.R,
          ctlL= x.L[idx],
          ctlIndL= x.indL[idx],
          ctlR= x.R[idx],
          ctlIndR= x.indR[idx],
          dist= dist[idx])
      }, .EACHI, on= c("breakL", "breakR"), nomatch= NULL]
      # }, .EACHI, on= c("breakL", "breakR"), nomatch= NULL]
      # for each (potentially duplicated) control pair, select closest enr pair
      setorderv(.c, "dist")
      .c <- .c[, .SD[1], .(ctlL, ctlR)]
      # Add control pairs to enr
      enr[.c, c("ctlL", "ctlIndL", "ctlR", "ctlIndR"):= .(i.ctlL, i.ctlIndL, i.ctlR, i.ctlIndR), on= c("L", "R")]
      # Remove control pairs from potential controls
      ctls[.c, used:= T, on= c("L==ctlL", "R==ctlR")]
      ctls <- ctls[is.na(used), !"used"]
      # Increase max dist when no more control pairs can be found
      if(!nrow(.c))
      {
        breaks <- breaks+0.5
        enr[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
        ctls[, c("breakL", "breakR"):= .(round(indL/breaks), round(indR/breaks))]
      }
      print(paste0(sum(!is.na(enr$ctlL)), "/", tot))
    }
    enr
  }, .(variable, min_counts)]
  res[dat, ctlLog2FoldChange:= i.log2FoldChange, on= c("ctlL==L", "ctlR==R")]
  res[vl_Dmel_motifs_DB_full, name:= i.motif_cluster, on="variable==motif_ID"]
  saveRDS(res,
          "db/activity_matched_controls/activity_matched_controls.rds")
}else
  res <- readRDS("db/activity_matched_controls/activity_matched_controls.rds")

# PLOT ----
pdf("pdf/draft/activity_predicitve_motifs_activity_matched.pdf", 
    height= 3.5,
    width= 7.5)
par(mar= c(2.1, 7, 2.1, 18),
    mgp= c(1.5,0.5,0),
    oma= c(0,0,0,5),
    las= 1,
    tcl= -0.2,
    xpd= NA)
res[actL!="Inactive" & actL!="Inactive", {
  xpos <- seq(1, 5, length.out= 3)
  vl_boxplot(ctlIndL,
             indL,
             ctlIndR,
             indR,
             ctlLog2FoldChange,
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
                   paste0("Motif pairs (>= ", min_counts, ")")),
         bty= "n")
  .SD
}, .(name, min_counts)]
dev.off()

file.show("pdf/draft/activity_predicitve_motifs_activity_matched.pdf")
