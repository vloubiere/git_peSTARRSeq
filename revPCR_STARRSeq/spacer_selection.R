setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(Biostrings)

# find NotI restriction sites gw 
dat <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)
dat <- dat[names(dat) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")]

dat <- vmatchPattern("GCGGCCGC", seq)
dat <- rbindlist(lapply(dat, as.data.table), idcol = T)
colnames(dat)[1] <- "seqnames"

# min and max coor
dat[, size := width(seq)[match(seqnames, names(seq))], seqnames]
dat[, c("min", "max") := .(round(size/10), round(size-(size/10))), seqnames]

# Overlap REs
peaks <- readRDS("Rdata/intron_STARRSeq/regulatory_elements_S2.rds")
dat[, no_RE:=  peaks[dat, all(abs(i.end-start)>3000 & i.start>min & i.end<max), .EACHI, on= "seqnames"]$V1]

# Overlap exons
TxDb <- "TxDb.Dmelanogaster.UCSC.dm3.ensGene"
.e <- as.data.table(exons(TxDb.Dmelanogaster.UCSC.dm3.ensGene, columns= "GENEID"))
.e <- .e[, .(GENEID= unlist(GENEID)), seqnames:strand]
.e <- .e[order(seqnames, start, end, GENEID)]
dat[, no_junction:= .e[dat, all(abs(i.start-start)>3000 & abs(i.start-end)>3000), .EACHI, on= "seqnames"]$V1]

# Resize
dat[, start:=  start-3000]
dat[, end:=  start+6000]

# Extract sequence
dat[(no_RE), seq:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(.SD)))]
dat[!is.na(seq), CG_perc:= length(which(grepl("C|G", unlist(strsplit(seq, "")))))/nchar(seq)*100, seq]
dat[, CG_check:= CG_perc>45 & CG_perc<55]

#-------------------------------------------------#
# All introns with good characteristics
#-------------------------------------------------#
regions <- dat[(no_RE) & (CG_check) & (no_junction)]

pdf("pdf/revPCR_STARRSeq/screenshot_potential_spacers.pdf", width = 20, height = 5)
par(mar= c(5,20,2,2))
regions[, {
  .c <- GRanges(.SD)
  .c <- unlist(GRangesList(.c, resize(.c, 20000, "center"), resize(.c, 100000, "center")))
  my_screenshot(bw_GR_list = c("../available_data/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                               "../available_data/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                               "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), names.arg= c("RNA-Seq", "ATAC-Seq", "DSCP_200bp", "DSCP_600bp", "RPS12_200bp", "RPS12_600bp"),
                bed = .c, gene.n = 20, enr_cutoff = c(50, 10, 25, 100, 100, 100))
}, seq(nrow(regions))]
dev.off()

final <- regions[c(2,3,4)]

pdf("pdf/revPCR_STARRSeq/screenshot_selected_spacers.pdf", width = 20, height = 5)
par(mar= c(5,20,2,2))
final[, {
  .c <- GRanges(.SD)
  .c <- unlist(GRangesList(.c, resize(.c, 20000, "center"), resize(.c, 100000, "center")))
  my_screenshot(bw_GR_list = c("../available_data/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                               "../available_data/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                               "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), names.arg= c("RNA-Seq", "ATAC-Seq", "DSCP_200bp", "DSCP_600bp", "RPS12_200bp", "RPS12_600bp"),
                bed = .c, gene.n = 20, enr_cutoff = c(50, 10, 25, 100, 100, 100))
}, seq(nrow(final))]
dev.off()
fwrite(final, "Rdata/revPCR_STARRSeq/spacers_sequence.txt", col.names = T, row.names = F, sep= "\t", quote= F)


