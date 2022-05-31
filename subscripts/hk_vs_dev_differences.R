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
# enhancers <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
hk <- fread("../gw_STARRSeq_bernardo/db/peaks/RPS12_600bp_gw_cut_merged.peaks.txt", 
            col.names = c("seqnames", "start"), 
            select = c(1,2))
# hk <- hk[V8>=3,. (seqnames= V1, start= V2, end= V2, enr= V7)]
dev <- fread("../gw_STARRSeq_bernardo/db/peaks/DSCP_600bp_gw_cut_merged.peaks.txt",
             col.names = c("seqnames", "start"), 
             select = c(1,2))
# dev <- dev[V8>=5,. (seqnames= V1, start= V2, end= V2, enr= V7)]
enhancers <- rbindlist(list(hk= hk, dev= dev), 
                       idcol = T)
enhancers[, end:= start]

# overlaps
dat <- vl_closestBed(enhancers, 
                     genes)
res <- dat[, .SD[which.min(abs(dist))], idx.a]
pl <- res[, .N, .(.id.a, FBgn.b)][, .N, .(.id.a, 
                                          cut(N, 
                                              c(0,1,2,5,Inf), 
                                              include.lowest= T,
                                              labels= c("1","2","3-5", ">5")))]
mat <- as.matrix(dcast(pl, cut~.id.a, value.var = "N"), 1)
norm <- apply(mat, 2, function(x) x/sum(x, na.rm= T)*100)
bar <- barplot(norm,
               col= RColorBrewer::brewer.pal(n = nrow(mat), 
                                             name = "Greens"),
               las= 2)

text(bar, 
     par("usr")[4],
     colSums(mat),
     xpd= T, 
     pos= 3)


test <- vl_closestBed(enhancers, 
                      enhancers, 
                      min_dist = 1)
test[, dist:= abs(dist)]
setorderv(test, "dist")
test <- test[.id.a==.id.b & is.finite(dist)]
test <- test[, mean(dist[1:5]), .(.id.a, idx.a)]
vl_boxplot(V1~.id.a, test, notch= T)

# GO <- unique(res[, .(.id.a, FBgn.b)])
# par(las= 1)
# pl <- vl_GO_clusters(split(GO$FBgn.b, GO$.id.a), plot= F)
# plot(pl, padj_cutoff = 1e-10)




