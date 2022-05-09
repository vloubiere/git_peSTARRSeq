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
if(!exists("RNASeq"))
{
  RNASeq <- data.table(file= list.files("db/FPKM/", full.names = T))
  RNASeq <- RNASeq[, fread(file), file]
  RNASeq <- RNASeq[, .(RNASeq_log2FPKM= log2(mean(V4)+0.001)), .(FBgn= V1)]
}

# Enhancer assignment
if(!exists("assignment"))
{
  enh <- unique(gtf[type=="gene" 
                    & seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"), 
                    .(seqnames, start, end, strand, FBgn= gene_id, symbol= gene_name)])
  enh[, enh_window_start:= ifelse(strand=="-", 
                                  end-2000,
                                  start-2000)]
  enh[strand=="-" & enh_window_start>end, enh_window_start:= end]
  enh[, enh_window_end:= ifelse(strand=="-", 
                                end+2000,
                                start+2000)]
  enh[strand=="+" & enh_window_end>end, enh_window_end:= end]
  setorderv(enh, c("seqnames", "start", "end", "FBgn"))
  enh[, gene_coor:= paste0(seqnames, ":", start, "-", end, ":", strand)]
  enh <- enh[, .(gene_coor, 
                 FBgn, 
                 symbol, 
                 seqnames, 
                 start= enh_window_start, 
                 end= enh_window_end)]
  assignment <- vl_closestBed(enh, STARR[, .(seqnames, start, end, ID, 
                                             hk_log2FoldChange, 
                                             hk_padj, 
                                             dev_log2FoldChange, 
                                             dev_padj, 
                                             enhancer_group)])
  assignment <- assignment[dist==0]
  setnames(assignment, new= function(x) gsub(".a$|.b$", "", x))
}

# Format data
dat <- assignment[hk_padj<0.01 & hk_log2FoldChange>1 | dev_padj<0.01 & dev_log2FoldChange>1, N:= .N, FBgn]
fwrite(fread("../gw_STARRSeq_bernardo/db/peaks/DSCP_200bp_gw.UMI_cut_merged.peaks.txt")[V8>5, .(seqnames=V1, start= V2-250, end= V2+250)],
       "db/peaks/DSCP_peaks_200bp_STARRSeq_Almeida.bed")
fwrite(fread("../gw_STARRSeq_bernardo/db/peaks/RpS12_200bp_gw.UMI_cut_merged.peaks.txt")[V8>5, .(seqnames=V1, start= V2-250, end= V2+250)],
       "db/peaks/RpS12_peaks_200bp_STARRSeq_Almeida.bed")
genes <- GRanges(dat[N==2, gene_coor])
genes <- unique(vl_resizeBed(genes, center = "center", upstream = end(genes)-start(genes), downstream = end(genes)-start(genes)))

vl_screenshot(bed = genes[6:10],
              tracks = c("../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                         "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw", 
                         "Rdata/BA_300bp_TWIST_STARRSeq.txt"), 
              genome = "dm3")



assignment[, bin:= floor(rleid(FBgn)/100)]
assignment[, length(unique(FBgn)), bin]
assignment[RNASeq, RNASeq_log2FPKM:= RNASeq_log2FPKM, on= "FBgn"]

number <- assignment[, .(N= length(unique(ID))), .(bin, FBgn, RNASeq_log2FPKM)]
plot(number[, .(median(N), median(RNASeq_log2FPKM, na.rm = T)), bin][, .(V1, V2)])






