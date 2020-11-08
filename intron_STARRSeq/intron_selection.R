require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(org.Dm.eg.db)
require(rtracklayer)

# genes
TxDb <- "TxDb.Dmelanogaster.UCSC.dm3.ensGene"
# transcripts
.e <- as.data.table(exons(TxDb.Dmelanogaster.UCSC.dm3.ensGene, columns= "TXNAME"))
.e <- .e[, .(TXNAME= unlist(TXNAME)), seqnames:strand]
.e <- .e[order(seqnames, start, end)]
.e <- na.omit(.e[, .(seqnames, start= end[-(.N)], end= start[-1], strand), .(seqnames, TXNAME, strand)])

.i <- unique(.e[, .(seqnames, start, end, strand)])

peaks <- data.table(file= list.files("../gw_STARRSeq_bernardo/db/peaks/", full.names = T))
peaks <- peaks[, fread(file, col.names = c("seqnames", "coor"), select = c(1,2)), file]

sel <- peaks[.i, !any(i.start<coor & i.end>coor), .EACHI, on= "seqnames"]$V1
sel <- .i[(sel)]
