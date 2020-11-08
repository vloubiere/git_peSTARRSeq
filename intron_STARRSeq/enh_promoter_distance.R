# setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
setwd("/home/vloubiere/projects/pe_STARRSeq/")
# sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
sapply(list.files("/home/vloubiere/functions/", ".R$", full.names = T), source)
options(datatable.print.topn= 3)
require(data.table)
require(rtracklayer)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm3)

if(!file.exists("Rdata/intron_STARRSeq/regulatory_elements_S2.rds"))
{
  #------------------------#
  # Compute TSSs active genes
  #------------------------#
  tss <- as.data.table(import("../../genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf"))
  tss <- tss[tss$type== "transcript" & seqnames %in% c("2L", "2R", "3L", "3R", "4", "X"), 
             .(seqnames= paste0("chr", seqnames), start= ifelse(strand=="+", start, end), end= ifelse(strand=="+", start, end), strand, gene_id, gene_name)]
  # S2 RNA-seq (from Arnold et al 2013
  .S2 <- fread("/groups/stark/gerlach/work/rnaseq_screen/data/RNAseq_S2_cut13-37_gene_rpkm.txt", fill= T, skip = 5)
  colnames(.S2) <- gsub("\\#", "", colnames(.S2))
  tss[.S2, S2_RNA:= i.V4, on= "gene_id==V1"]
  tss <- unique(tss)
  
  #------------------------#
  # Compute ATAC-Seq peak
  #------------------------#
  binsize <- 125
  bins <- as.data.table(getChromInfoFromUCSC("dm3"))
  bins <- bins[c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"), , on= "chrom"]
  bins <- bins[, .(seqnames= chrom, start= seq(1000, length, binsize)), (bins)]
  bins[, end:= start+binsize-1]
  bins <- bins[, .(seqnames, start, end)]
  .q <- data.table(file= list.files("../available_data/db/bed/", "GSE119708_ATAC", full.names = T))
  .q <- .q[, my_countReads(bins, file, sorted = T), file]
  .q[, enr:= log2(counts+1)-log2(mean(.SD$counts, na.rm= T)+1), file]
  sel <- .q[, rep(all(enr>1.5), .N), .(seqnames, start, end)]$V1
  # Merge peaks
  ATAC <- GRanges(unique(.q[(sel), seqnames:end]))
  ATAC <- as.data.table(reduce(ATAC, min.gapwidth= 2*binsize))
  
  #---------------------------#
  # STARR-Seq 
  #---------------------------#
  STARR <- data.table(file= list.files("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/peaks/", full.names = T))
  STARR[, type:= gsub("RpS", "RPS", tstrsplit(basename(file), "bp_gw", keep = 1)), file]
  STARR <- STARR[, fread(file, select = c(1,2,8,9), col.names = c("seqnames", "start", "enr", "qvalue")), (STARR)]
  STARR[, end:= start]
  
  #---------------------------#
  # All peaks
  #---------------------------#
  tss$type <- "tss"
  ATAC$type <- "ATAC"
  peaks <- rbind(tss, ATAC, STARR, fill= T)
  peaks <- peaks[, .(seqnames, start, end, strand= ifelse(is.na(strand), "*", as.character(strand)), type, gene_id, gene_name, S2_RNA, enr, qvalue)] 
  # Final filtering
  peaks <- peaks[type=="tss" & S2_RNA>=2 | type!="tss"]
  peaks[, c("start", "end"):= .(round(rowMeans(.SD)), round(rowMeans(.SD))), .SDcols= c("start", "end")]
  saveRDS(peaks, "Rdata/intron_STARRSeq/regulatory_elements_S2.rds")
}
peaks <- readRDS("Rdata/intron_STARRSeq/regulatory_elements_S2.rds")

# Measure distance to closest
peaks$c_tss <- peaks[type=="tss"][peaks, min(abs(i.start-start[!(type==i.type & start==i.start)])), .EACHI, on= "seqnames"]$V1

pdf("pdf/intron_STARRSeq/enh_promoter_distance.pdf", width = 8)
par(mar= c(5,8,2,2))
my_boxplot(c_tss~type, peaks, bw= 1000, horizontal = T, ylim= c(0, 25e3), las= 1, outline= T, xaxt= "n", xlab= "closest active TSS")
axis(1, c(0, 2e3, 5e3, 10e3, 25e3))
abline(v=2e3, col= "red")
dev.off()





