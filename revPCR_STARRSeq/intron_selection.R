setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
# setwd("/home/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
# sapply(list.files("/home/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(org.Dm.eg.db)
require(rtracklayer)

if(!exists("sel"))
{
  # genes
  TxDb <- "TxDb.Dmelanogaster.UCSC.dm3.ensGene"
  # exons
  .e <- as.data.table(exons(TxDb.Dmelanogaster.UCSC.dm3.ensGene, columns= "GENEID"))
  .e <- .e[, .(GENEID= unlist(GENEID)), seqnames:strand]
  .e <- .e[order(seqnames, start, end, GENEID)]
  # introns
  .i <- na.omit(.e[, .(seqnames, start= end[-(.N)], end= start[-1], strand), .(seqnames, GENEID, strand)])
  .i <- unique(.i[, .(seqnames, start, end, strand, GENEID)])
  # Compute overlaps
  .i$intron <- .i[.i, any(between(start, i.start, i.end, incbounds = F)|between(end, i.start, i.end, incbounds = F)), .EACHI, on= "seqnames"]$V1
  peaks <- readRDS("Rdata/intron_STARRSeq/regulatory_elements_S2.rds")
  .i$peaks <- peaks[.i, any(between(start, i.start, i.end, incbounds = F)|between(end, i.start, i.end, incbounds = F)), .EACHI, on= "seqnames"]$V1
  .i[peaks, expression:= S2_RNA,  on= "GENEID==gene_id"]
  # Selection 
  sel <- .i[end-start > 1750 & end-start < 2250 & expression>20 & !intron & !peaks]
}

#-------------------------------------------------#
# All introns with good characteristics
#-------------------------------------------------#
regions <- unique(sel[, .(seqnames, start, end)])
regions <- regions[, .(seqnames, start= c(start-2e3, start-10e3), end= c(end+2e3, end+10e3), name= .GRP), regions][, 4:7]

pdf("pdf/intron_STARRSeq/screenshot_introns.pdf", width = 15)
par(mar= c(5,20,2,2))
regions[, {
  .c <- GRanges(.SD)
  my_screenshot(bw_GR_list = c("../available_data/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                               "../available_data/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                       "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                       "../gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                       "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                       "../gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), bed = .c, gene.n = 20, enr_cutoff = c(100, 20, 50, 150, 150, 150))
}, name]
dev.off()

#-------------------------------------------------#
# Selected introns
#-------------------------------------------------#
sub <- regions[name %in% c(2,6,23,47,52,71)]              
pdf("pdf/intron_STARRSeq/screenshot_selected_introns.pdf", width = 15, height = 5)
par(mar= c(5,20,2,2))
sub[, {
  .c <- GRanges(.SD)
  my_screenshot(bw_GR_list = c("../available_data/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                               "../available_data/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                               "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                               "../gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), names.arg= c("RNA-Seq", "ATAC-Seq", "DSCP_200bp", "DSCP_600bp", "RPS12_200bp", "RPS12_600bp"),
                bed = .c, gene.n = 20, enr_cutoff = c(100, 20, 50, 150, 150, 150))
}, name]
dev.off()

#-------------------------------------------------#
# Schnurri locus
#-------------------------------------------------#
pdf("pdf/intron_STARRSeq/screenshot_schnurri.pdf", width = 10, height = 3)
par(mar= c(2,10,2,2))
my_screenshot(bw_GR_list = c("../available_data/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                             "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw"), 
              bed = GRanges("chr2R", IRanges(7048428, 7111206)), names.arg= c("RNA-Seq", "DSCP_200bp"),
              gene.n = 20, enr_cutoff = c(30, 100))
dev.off()





