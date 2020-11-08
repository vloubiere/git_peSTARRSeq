# setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
setwd("/home/vloubiere/projects/pe_STARRSeq/")
# sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
sapply(list.files("/home/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(org.Dm.eg.db)
require(rtracklayer)

# genes
TxDb <- "TxDb.Dmelanogaster.UCSC.dm3.ensGene"
# exons
.e <- as.data.table(exons(TxDb.Dmelanogaster.UCSC.dm3.ensGene, columns= "TXNAME"))
.e <- .e[, .(TXNAME= unlist(TXNAME)), seqnames:strand]
.e <- .e[order(seqnames, start, end)]
# introns
.i <- na.omit(.e[, .(seqnames, start= end[-(.N)], end= start[-1], strand), .(seqnames, TXNAME, strand)])
.i <- unique(.i[, .(seqnames, start, end, strand)])

# Check

if(!exists(".g"))
{
  .g <- as.data.table(import("../../genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf"))
}
# region <- resize(GRanges(.g[type=="gene" & gene_name=="RpS12"]), 10000, "center")
# region <- resize(GRanges(.g[type=="gene" & gene_name=="N"]), 100000, "center")
region <- resize(GRanges(.g[type=="gene" & gene_name=="shn"]), 50000, "center")
seqlevelsStyle(region) <- "UCSC"
.p <- GRangesList(GRanges(NULL), GRanges(NULL), GRanges(NULL), GRanges(NULL), GRanges(NULL), 
                  GRanges(.e[, seqnames:end]), 
                  GRanges(.i[, seqnames:end]))

pdf("pdf/intron_STARRSeq/screenshot_selected_intron.pdf", width = 15)
par(mar= c(5,20,2,2))
my_screenshot(bw = c("../../public/RNASeq_S2/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                     "../../public/PROSeq_s2/GSM2055132_S2_NHS_PRO_seq_norm_plus.bw",
                     "../../public/PROSeq_s2/GSM2055132_S2_NHS_PRO_seq_norm_minus.bw",
                     "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                     "../gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                     "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                     "../gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), bed = region, peaks = .p)
dev.off()
              
