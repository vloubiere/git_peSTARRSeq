setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)


# Import genes
genes <- rtracklayer::import("../../genomes/Drosophila_melanogaster/flybase/dm3/dmel-all-filtered-r5.57_genes_and_transcripts_only.gff")
GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[Name %in% c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo")]
genes <- vl_resizeBed(genes,
                      "start",
                      upstream = 10000,
                      downstream = genes[,end-start+1])

# Import lib and overlap
lib <- readRDS("Rdata/vl_library_twist008_112019.rds")
lib <- as.data.table(lib)
lib[genes, gene:= i.Name, on= c("seqnames", "end>=start", "start<=end")]

# Import dat and add left an right genes
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
dat[lib, geneL:= i.gene, on="L==ID_vl"]
dat[lib, geneR:= i.gene, on="R==ID_vl"]
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Find best activity-matched controls ----
dat$ctlL <- dat$ctlR <- NULL
enr <- dat[geneL==geneR]
ctls <- dat[geneL!=geneR | is.na(geneL) | is.na(geneR)]
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
enr[dat, ctlLog2FoldChange:= i.log2FoldChange, on= c("ctlL==L", "ctlR==R")]

pdf("pdf/draft/multi_enhancer_genes_vs_controls.pdf", 3.5, 3)
par(tcl= -0.2,
    las= 1,
    mgp= c(2,.5,0),
    mar= c(5,4,2,8))
xpos <- seq(1, 5, length.out= 3)
vl_boxplot(enr[, .(ctlIndL, indL, ctlIndR, indR, ctlLog2FoldChange, log2FoldChange)],
           compute_pval= list(c(1,2), c(3,4), c(5,6)),
           notch= T,
           xaxt= "n",
           col= c("lightgrey", "rosybrown1"),
           ylab= "Activity (log2)",
           at= rep(xpos, each= 2)+c(-0.35, 0.35))
axis(1, at= xpos, labels = c("5'", "3'", "Pair"))
legend(par("usr")[2],
       par("usr")[4],
       fill= c("lightgrey", "rosybrown1"),
       legend= c("Controls (diff. loci)",
                 "Same locus"),
       bty= "n",
       xpd= NA,
       cex= 0.7)
vl_boxplot(enr[actL!="Inactive" & actR!="Inactive",
               .(ctlIndL, indL, ctlIndR, indR, ctlLog2FoldChange, log2FoldChange)],
           compute_pval= list(c(1,2), c(3,4), c(5,6)),
           notch= T,
           xaxt= "n",
           col= c("lightgrey", "rosybrown1"),
           ylab= "Activity (log2)",
           at= rep(xpos, each= 2)+c(-0.35, 0.35),
           main= "Active pairs only")
axis(1, at= xpos, labels = c("5'", "3'", "Pair"))
legend(par("usr")[2],
       par("usr")[4],
       fill= c("lightgrey", "rosybrown1"),
       legend= c("Controls (diff. loci)",
                 "Same locus"),
       bty= "n",
       xpd= NA,
       cex= 0.7)
dev.off()