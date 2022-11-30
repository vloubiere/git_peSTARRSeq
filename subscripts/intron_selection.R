setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

#-------------------------------#
# Make intron objects
#-------------------------------#
if(!file.exists("/groups/stark/vloubiere/genomes/flybase/dm3/dm3_introns.txt"))
{
  dat <- as.data.table(exons(TxDb.Dmelanogaster.UCSC.dm3.ensGene, columns= c("TXNAME")))
  dat <- dat[, .(FBtr= unlist(TXNAME)), setdiff(colnames(dat), "TXNAME")]
  setorderv(dat, c("FBtr", "start"))
  introns <- dat[, .(seqnames= rep(seqnames[1], .N-1), start= end[-(.N)]+1, end= start[-1]-1, strand= strand[1]), FBtr]
  # ADD FBgn
  .g <- as.data.table(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene, columns= c("TXNAME", "GENEID")))
  .g[, GENEID:= unlist(GENEID)]
  .g <- .g[, .(FBtr= unlist(TXNAME)), setdiff(colnames(.g), "TXNAME")]
  introns[.g, FBgn:= i.GENEID, on= "FBtr"]
  fwrite(introns, 
         "/groups/stark/vloubiere/genomes/flybase/dm3/dm3_introns.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}else
  introns <- fread("/groups/stark/vloubiere/genomes/flybase/dm3/dm3_introns.txt", 
                   colClasses = c("character", "character", "numeric", "numeric", "character", "character"))

vl_screenshot(data.table(seqnames= "chrX", start= 3026905-20000, end= 3068295+20000), 
              "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE119708_ATAC_rep1_uniq.bw", 
              genome= "dm3",
              highlight_regions = introns)

#-----------------------#
# Selection 
#-----------------------#
if(!file.exists("db/library_design/alternative_spacers/intron_candidates.txt"))
{
  # 1- No overlap with other junctions
  sel <- unique(introns[, .(seqnames, start, end, strand, FBgn)])
  all_junct <- sel[, .(seqnames, start= unlist(.(start_old, end_old))), .(seqnames, start_old= start, end_old= end)]
  all_junct <- unique(all_junct[, .(seqnames, start, end= start)])
  sel <- sel[(all_junct[sel, .N==0, .EACHI, on= c("seqnames", "start<end", "end>start")]$V1)]
  # 2- No overlap with any predicted promoter
  prom <- as.data.table(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene))
  prom[strand== "+", c("start", "end"):= .(start-1000, start+1000)]
  prom[strand== "-", c("start", "end"):= .(end-1000, end+1000)]
  sel <- sel[(prom[sel, .N==0, .EACHI, on= c("seqnames", "start<end", "end>start")]$V1)]
  # 3- SIZE
  sel <- sel[between(end-start, 275, 325) 
             | between(end-start, 1750, 2250) 
             | between(end-start, 4750, 5250) 
             | between(end-start, 8000, 12000)]
  # 4- EXPRESSION
  exp <- fread('../../genomes/flybase/dm6/gene_rpkm_report_fb_2017_05.tsv', fill= T, skip= 5)
  sel <- sel[FBgn %in% exp[RNASource_name=="mE_mRNA_S2R+_cells" & Bin_value>1, `FBgn#`]]
  # 5- Does not overlap with any peak (ATAC, DHS, ATAC-Seq, H3K4me3)
  peaks <- data.table(file= c("db/peaks/ATAC_peaks.txt",
                              "db/peaks/DHS_peaks.txt",
                              "db/peaks/H3K4me3_peaks.txt", 
                              "db/peaks/STARR_DSCP_200_peaks.txt",
                              "db/peaks/STARR_DSCP_600_peaks.txt",
                              "db/peaks/STARR_RPS12_200_peaks.txt",
                              "db/peaks/STARR_RPS12_600_peaks.txt"))
  peaks <- peaks[, fread(file, colClasses = c("character", rep("numeric", 5))), file]
  sel <- sel[(peaks[sel, .N==0, .EACHI, on= c("seqnames", "start<end", "end>start")]$V1)]
  
  #SAVE
  sel <- sel[order(end-start)]
  sel[end-start<=325, size:= 300]
  sel[between(end-start, 1750, 2250), size:= 2000]
  sel[between(end-start, 4750, 5250), size:= 5000]
  sel[between(end-start, 8000, 12000), size:= 10000]
  fwrite(sel, "db/library_design/alternative_spacers/intron_candidates.txt", col.names = T, na= NA)
}else
  sel <- fread("db/library_design/alternative_spacers/intron_candidates.txt")

#-------------------------------------------------#
# Screenshots
#-------------------------------------------------#
if(!file.exists("pdf/design/screenshot_all_intron_candidates.pdf"))
{
  sel$sub <- rep(seq(100), each= 5)[1:nrow(sel)]
  
  pdf("pdf/design/screenshot_all_intron_candidates.pdf", width = 15)
  sel[, {
    .c <- GRanges(.SD)
    vl_screenshot(c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                    "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                    "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                    "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                    "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                    "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), 
                  genome= "dm3",
                  bed = resize(.c, end-start+(2*size), "center") , 
                  max = c(100, 20, 50, 150, 150, 150), 
                  highlight_regions = .c)
    mtext(paste(size, "sel=", paste0(.I, collapse= "-")))
  }, .(size, sub)]
  dev.off()
}

#-------------------------------------------------#
# Selected introns
#-------------------------------------------------#
sub <- sel[c(3,15,48,257,264,291,348,350,357,365,380,381), !"sub"]              

pdf("pdf/design/screenshot_selected_introns.pdf", width = 15)
sub[, {
  .c <- GRanges(seqnames = seqnames, IRanges(start, end))
  vl_screenshot(c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSM480160_GA0840_Drosophila_S2_RNAseq.bw",
                  "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE119708_ATAC_rep1_uniq.bw",
                  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
                  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw",
                  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/RPS12_600bp_gw_cut_merged.bw"), 
                genome= "dm3",
                bed = resize(.c, end-start+(2*size), "center") ,
                max = c(100, 20, 50, 150, 150, 150), 
                highlight_regions = .c)
  mtext(paste(size, "sel=", paste0(.I, collapse= "-")))
}, (sub)]
dev.off()

fwrite(sub, 
       "db/library_design/alternative_spacers/selected_introns.txt")
