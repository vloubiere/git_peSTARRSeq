setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(GenomicRanges)

# Import DHS Data
DHS <- fread("/groups/stark/almeida/Projects/Sequence_rules_epigenomic_layers/results/20220125_Enh_act_and_chromatin/DHS_peaks_annotated_enh_act.txt")
DHS <- DHS[Enh_short_overlap=="Inactive"]

# Import STARR-Seq peaks
load_peak_function <- function(x, width){
  path<- "/groups/stark/almeida/data/STARRseq/dm3/20190401_200bp_gw_STARRseq/data/"
  # load peaks file and convert to granges
  gr <- makeGRangesFromDataFrame(read.delim(paste0(path, x), header = F),
                                 seqnames.field = "V1", start.field = "V2", end.field = "V2", keep.extra.columns = T)
  # resize peaks
  gr <- resize(gr, width, "center")
  # Treat peak information
  names(mcols(gr))[5:7] <- c("Enrch.", "Corr_enrch", "p_value")
  as.data.table(gr)
}
dCP_200bp <- load_peak_function("DSCP_200bp_gw.UMI_cut_merged.peaks.txt", 249)
dCP_200bp <- dCP_200bp[Corr_enrch>3]

# Select DHS peaks that do not overlap dev STARR-Seq or promoter
DHS <- DHS[dCP_200bp[DHS, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N==0]
DHS[, summit:= start+summit]
DHS[, start:= summit-124]
DHS[, end:= summit+124]
prom <- rtracklayer::import("../../genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf")
seqlevelsStyle(prom) <- "UCSC"
prom <- as.data.table(prom)
prom <- vl_resizeBed(prom[type=="transcript"], "start", 100, 50)
DHS <- DHS[prom[DHS, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N==0]
setorderv(DHS, "fold_enrichment")

# DHS enrichment cutoff
DHS_tmp <- tempfile(fileext = ".bed")
rtracklayer::export(DHS, DHS_tmp)
STARR_tmp <- tempfile(fileext = ".bed")
rtracklayer::export(dCP_200bp, STARR_tmp)
pl <- as.call(list(quote(vl_screenshot),
                   bed= NULL,
                   tracks = c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw", STARR_tmp,
                              "/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/tracks/S2_DHS.bw", DHS_tmp),
                   names = c("STARR-Seq", "STARR-Seq peaks", "DHS", "DHS peaks")))
cuts <- seq(min(DHS$fold_enrichment), max(DHS$fold_enrichment), length.out= 10)
idx <- sapply(cuts, function(x) which.min(abs(DHS[-c(1,2), fold_enrichment]-x)))

par(mfrow=c(3,1))
pl[[2]] <- DHS[idx, .(seqnames, start= start-2500, end= end+2500)]
eval(pl)
pl[[2]] <- DHS[idx+1, .(seqnames, start= start-2500, end= end+2500)]
eval(pl)
pl[[2]] <- DHS[idx+2, .(seqnames, start= start-2500, end= end+2500)]
eval(pl)

DHS <- DHS[fold_enrichment>=3] 

# Select rolling window with max number of DHS sites and twist motifs
dat <- rbind(dCP_200bp[, .(seqnames, start, end, strand, type= "STARR")],
             DHS[, .(seqnames, start, end, strand, type= "DHS")])
dat <- dat[seqnames %in% c("chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4")]
setorderv(dat, c("seqnames", "start", "end"))
counts <- vl_motif_counts(dat, 
                          sel= "cisbp__M2014",
                          genome= "dm3", 
                          collapse_overlapping= T,
                          p.cutoff = 1e-4)
dat$twist_counts <- counts[[1]]
dat[, idx:= rowid(seqnames)]
dat[, sel_start:= idx-425]
dat[, sel_end:= idx+424]
dat$N_DHS <- dat[dat, sum(type=="DHS"), .EACHI, on= c("seqnames", "idx<=sel_end", "idx>=sel_start")]$V1
dat$N_twist <- dat[dat, sum(between(twist_counts, 2, 5)), .EACHI, on= c("seqnames", "idx<=sel_end", "idx>=sel_start")]$V1

plot(dat[order(N_twist), N_twist], type= "l")
lines(dat[order(N_twist), N_DHS], type= "l", col= "red")
legend("topleft", 
       col= c("black", "red"),
       lty= 1,
       legend= c("twist", "DHS"))
abline(h= 63, col= "red")
abline(h= 62, col= "black")

dat[N_DHS>61 & N_twist>65]
dat[, sel:= {
  idx <- which(seqnames=="chr2R"  & start==9183768 & end==9184016)
  res <- rep(FALSE, .N)
  if(length(idx==1))
    res[(idx-425):(idx+424)] <- TRUE
  res
}, seqnames]
sel <- dat[(sel)]
table(sel[, type])
table(sel[, sum(twist_counts>=2)])

# Select control sequences from twist
TWIST <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
TWISTctls <- TWIST[hk_padj>0.05 & dev_padj>0.05 
              & dev_log2FoldChange<1 & hk_log2FoldChange<1 
              & seqnames==sel$seqnames[1] & start>=min(sel$start) & end<=max(sel$end)]
# Add random ones to reach 150 total
STARR_200bp <- rbind(load_peak_function("DSCP_200bp_gw.UMI_cut_merged.peaks.txt", 249),
                     load_peak_function("RpS12_200bp_gw.UMI_cut_merged.peaks.txt", 249))
set.seed(1)
rdm <- vl_random_regions_BSgenome("dm3", n = 2.5e5, width = 249, restrict_seqnames = sel$seqnames[1])
rdm <- rdm[vl_covBed(rdm, STARR_200bp)==0
           & vl_covBed(rdm, TWIST)==0
           & vl_covBed(rdm, DHS)==0
           & as.character(seqnames)==sel$seqnames[1] 
           & start>=min(sel$start)
           & end<=max(sel$end)]
set.seed(1)
rdm <- rdm[sample(nrow(rdm), 150-nrow(TWISTctls))]

#-----------------#
# Final
#-----------------#
lib <- rbindlist(list(sel[, .(type, seqnames, start, end)],
                      TWISTctls[, .(type= "control", seqnames, start, end)],
                 rdm[, .(type= "control", seqnames, start, end)]))
lib[, strand:= "+"]
lib[, enh_sequence:= vl_getSequence(.SD, genome = "dm3")]
saveRDS(lib, 
        "db/library_design/twist015/DHS_sublib_twist015.rds")

