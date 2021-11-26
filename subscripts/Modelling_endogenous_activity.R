setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

#------------------------------------------#
# Import Vanja TSSs with activity
#------------------------------------------#
tss <- read_xlsx("/groups/stark/vloubiere/projects/PROSeq_leo/Rdata/COF.STAPseq.TSSs.clustered.with.metadata.xlsx")
tss <- as.data.table(tss)
tss <- tss[S2_RNAseq_RPKM!="NA" & S2_whole_cell_CAGE!="NA"]
cols <- names(tss)[-c(1,2,6,7,8,9,10,11,33)]
tss[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]
tss <- tss[, .(FBgn= nearest.gene.sense.flybase_id, 
               FBtr= nearest.transcript.sense,
               chr,
               start, 
               end, 
               TSS,
               strand,
               TSS_cluster,
               RNASeq_log2_FPKM_Vanja= log2(S2_RNAseq_RPKM+min(S2_RNAseq_RPKM[S2_RNAseq_RPKM>0])),
               CAGE_log2_FPKM= log2(S2_whole_cell_CAGE+min(S2_whole_cell_CAGE[S2_whole_cell_CAGE>0])))]

#------------------------------------------#
# Import STAP-Seq
#------------------------------------------#
STAP <- readRDS("/groups/stark/lorbeer/vincent/STAP_oligo.RDS")
tss[, STAP_id:= paste0(chr, "_", start, "_", end, "_", strand)]
STAP[tss, FBgn:= i.FBgn, on=  "oligo_id==STAP_id"]
STAP <- STAP[, .(FBgn, 
                 STAP_zfh1_log2FC= log.zfh1_responsiveness, 
                 STAP_ssp3_log2FC= log2((ssp3_enh_Rep1+1)/sum(ssp3_enh_Rep1)))]

#------------------------------------------#
# Import PRO-Seq
#------------------------------------------#
annotation <- fread("/groups/stark/hendy/Projects/Chromatin_remodeler_specificity/GEO_submission/GEO_submission_20210902/PROseq.tss.annotation.tsv")
annotation[, tss:= paste0(seqnames, ":", start, "-", end, ":", strand)]
PROSeq_counts <- fread("/groups/stark/hendy/Projects/Chromatin_remodeler_specificity/GEO_submission/GEO_submission_20210914/PROseq.tss.counts.tsv")
PROSeq_counts[annotation, c("FBgn", "symbol"):= .(i.gene_id, i.gene_name), on= "tss"]
cols <- c("Parental_0hrIAA_1", "Parental_0hrIAA_2")
PROSeq_counts[, paste0(cols, "_norm_counts") := lapply(.SD, function(x) (x+1)/sum(x)*1e6), .SDcols= cols]
PROSeq_counts[, PROSeq_log2_TSS_rpm := log2(rowMeans(.SD)), .SDcols= patterns("norm_counts$")]
PROSeq_counts <- PROSeq_counts[, .SD[which.max(PROSeq_log2_TSS_rpm)], .(FBgn, symbol), 
                               .SDcols= c("tss", "PROSeq_log2_TSS_rpm")]

#------------------------------------------#
# Import FLYBASE RNA-Seq
#------------------------------------------#
fpkm <- fread("/groups/stark/vloubiere/genomes/flybase/dm6/gene_rpkm_report_fb_2017_05.tsv", 
              fill= T,
              skip = 5)
fpkm <- fpkm[RNASource_name=="mE_mRNA_S2R+_cells", .(FBgn= `FBgn#`, RNASeq_log2_RPKM_flybase= log2(RPKM_value+1))]

#------------------------------------------#
# Make data table
#------------------------------------------#
dat <- merge(tss,
             PROSeq_counts[, .(FBgn, symbol, tss, PROSeq_log2_TSS_rpm)], 
             by= "FBgn")
dat <- merge(dat, 
             STAP)
dat <- merge(dat,
             fpkm)
plot(dat$STAP_zfh1_log2FC, dat$STAP_ssp3_log2FC)
