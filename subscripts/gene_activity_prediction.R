setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Functions
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import data
if(!exists("gtf"))
{
  gtf <- rtracklayer::import("/groups/stark/haberle/general.data/gene.annotation/dm3/Drosophila_melanogaster.BDGP5.78.chr.names.gtf")
  gtf <- as.data.table(gtf)[type=="gene"]
}
if(!exists("basal"))
{
  gene_correspondance <- loadRData("/groups/stark/haberle/general.data/gene.sets/dm3/Unique.annotated.TSSs.from.FlyBase.dmel-all-filtered-r5.57.corrected.by.S2.and.developmental.CAGE.in.windows.250.and.500bp.RData")
  gene_correspondance <- as.data.table(gene_correspondance)
  gene_correspondance <- gene_correspondance[, .SD[which.max(CAGE_tpm), .(seqnames, start= CAGE_TSS, end= CAGE_TSS, strand)], .(FBgn= gene_id, symbol= gene_name)]
  gene_correspondance <- vl_resizeBed(gene_correspondance, center = "start", upstream = 67, downstream = 66)
  gene_correspondance[, oligo_id:= paste0(seqnames, "_", start, "_", end, "_", strand)]
  basal <- fread("../STAPSeq_vanja/db/normalized_counts/GSE116197_Enhancer_screens_Normalized_tagcounts_per_oligo.txt.gz")
  basal <- merge(basal[, .(oligo_id, zfh1_enh_Rep1, zfh1_enh_Rep2, ssp3_enh_Rep1, no_enh_Rep1)],
                 gene_correspondance,
                 by= "oligo_id")
  basal <- basal[, .(FBgn, zfh1_enh_Rep1, zfh1_enh_Rep2, ssp3_enh_Rep1, no_enh_Rep1)]
}
if(!exists("PROSeq"))
{
  PROSeq <- lapply(list.files("/groups/stark/haberle/data/public/dm3/S2_cells/PROseq/GSE81649/processed/control.samples/", 
                              full.names = T), fread)
  PROSeq <- rbindlist(PROSeq, idcol = T)
  PROSeq[, .id:= paste0("PROSeq_rep", .id)]
  PROSeq <- dcast(PROSeq, V1~.id, value.var = "V2")
  setnames(PROSeq, "V1", "FBgn")
  cols <- names(PROSeq)[-1]
  PROSeq[, (cols):= lapply(.SD, function(x) (x+1)/sum(x+1)*1e6), .SDcols= cols]
  PROSeq[, mean_norm_counts:= rowMeans(.SD), .SDcols= patterns("^PROSeq")]
  PROSeq[gtf, gene_length:= end-start+1, on= "FBgn==gene_id"]
  PROSeq[, PROSeq_log2FPKM:= log2(mean_norm_counts/gene_length*1000)]
}
if(!exists("STARR-Seq"))
  STARR <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
if(!exists("dev"))
{
  dev <- fread("../gw_STARRSeq_bernardo/db/peaks/DSCP_200bp_gw.UMI_cut_merged.peaks.txt")
  dev <- dev[V8>=2,. (seqnames= V1, start= V2, end= V2, enr= V7)]
}
if(!exists("hk"))
{
  hk <- fread("../gw_STARRSeq_bernardo/db/peaks/RPS12_600bp_gw_cut_merged.peaks.txt")
  hk <- hk[V8>=1,. (seqnames= V1, start= V2, end= V2, enr= V7)]
} 
if(!exists("RNASeq"))
{
  RNASeq <- data.table(file= list.files("db/FPKM/", full.names = T))
  RNASeq <- RNASeq[, fread(file), file]
  RNASeq <- RNASeq[, .(RNASeq_log2FPKM= log2(mean(V4)+0.001)), .(FBgn= V1)]
}
if(!exists("genes"))
{
  genes <- unique(gtf[type=="gene" 
                      & seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"), 
                      .(seqnames, start, end, strand, FBgn= gene_id, symbol= gene_name)])
  setorderv(genes, c("seqnames", "start", "end"))
  # genes[, bin:= floor(rleid(FBgn)/100), seqnames]
  # genes[, bin:= .GRP, .(seqnames, bin)]
}
if(!exists("proms"))
  proms <- vl_resizeBed(genes, center = "start", upstream = 2000, downstream = 2000)
if(!exists("windows"))
  windows <- vl_resizeBed(genes, center = "start", upstream = 15000, downstream = genes[, end-start]+5000, ignore.strand = F)[]

# Quantif bins
dat <- windows[RNASeq, RNASeq_log2FPKM:= RNASeq_log2FPKM, on= "FBgn"]
dat <- cbind(dat,
             dev[dat, .(N_enh= .N, sum_enh= log2(sum(enr, na.rm= T))), .EACHI, on= c("seqnames", "start<=end", "end>=start")])
dat <- dat[N_enh>0 & !is.na(RNASeq_log2FPKM)]
dat <- dat[N_enh>0 & !is.na(RNASeq_log2FPKM)]
dat[, bin:= floor(rleid(FBgn)/100), seqnames]
dat[, bin:= .GRP, .(seqnames, bin)]
pls <- list(`Num. of enhancers (median)`= dat[, .(median(N_enh), median(RNASeq_log2FPKM)), bin][, !"bin"],
            `Sum of enhancer strenghts (median of bin) [log2]`= dat[, .(median(sum_enh), median(RNASeq_log2FPKM)), bin][, !"bin"])

par(mfrow= c(2,2),
    las= 1)
for(i in seq(pls))
{
  .c <- pls[[i]] 
  plot(.c,
       xlab= names(pls)[i],
       ylab= "RNA-Seq RPKM (log2)")
  legend("bottomright",
         legend= paste0("PCC=", round(cor.test(.c$V1, .c$V2)$estimate, 2)),
         bty= "n")
}
