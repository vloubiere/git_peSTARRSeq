setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# IMPORT GENES
if(!exists("proms"))
{
  genes <- fread("../../genomes/flybase/dm6/hk_dev_assignment_genes_fb_2017_05.txt")
  gtf <- rtracklayer::import("/groups/stark/annotations/dm3/dmel-all-filtered-r5.57_genes_and_transcripts_only.gff")
  seqlevelsStyle(gtf) <- "UCSC"
  gtf <- as.data.table(gtf)
  genes[gtf[type=="gene"], c("seqnames", "start", "end", "strand"):= .(i.seqnames, i.start, i.end, i.strand), on= "FBgn==ID"]
  genes <- na.omit(genes)
  genes[, length:= end-start+1]
  proms <- vl_resizeBed(genes, center = "start", upstream = 0, downstream = 0)
  # Keep only genes expressed in S2
  # expressed <- readRDS("/groups/stark/pleyer/projects/2021_CP_different_RNAPII_recruitment_strategies/202103_PROseq_complexII_mutatedBCs/analysis/differential_expression/expressedgenes/data/S2_expressedgenes.rds")
  # proms <- proms[FBgn %in% expressed$TSS$control]
}

# IMPORT ENHANCERS
dev <- fread("../gw_STARRSeq_bernardo/db/peaks/DSCP_200bp_gw.UMI_cut_merged.peaks.txt",
             col.names = c("seqnames", "start"), 
             select = c(1,2))
dev[, end:= start]

# overlaps
dat <- vl_closestBed(dev, 
                     dev, 
                     min_dist = 1)
enh_enh <- log10(dat[, .SD[which.min(abs(dist))], idx.a][, abs(dist)]+1)

pdf("pdf/draft/enh_enh_distance.pdf",
    height = 3,
    width = 2)
par(mar= c(3,5,1,2),
    mgp= c(3.5, 0.5, 0),
    tcl= -0.2)
vl_boxplot(list("devEnhancers"= enh_enh), 
           notch= T, 
           violin= T,
           yaxt= "n",
           ylim= c(1,5))
axis(2, 
     at= c(1,2,3,4, 5), 
     labels = c(10,100,"1,000", "10,000", "100,000"),
     las= 1)
title(ylab= "Distance to closest enhancer")
abline(h= log10(c(300, 2000)),
       lty= 2)
text(par("usr")[2],
     log10(c(300, 2000)),
     c("300bp", "2kb"),
     xpd= T,
     pos= 4,
     cex= 0.5)
dev.off()
