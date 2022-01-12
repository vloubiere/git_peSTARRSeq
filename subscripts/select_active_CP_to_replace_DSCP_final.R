setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(BSgenome.Dmelanogaster.UCSC.dm3)

#-------------------#
# Import 72k STAP-Seq library
#-------------------#
lib <- readRDS("/groups/stark/lorbeer/vincent/STAP_oligo.RDS")
lib[, sequence:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(seqnames, IRanges(start, end), strand)))]
# add DSCP
DSCP <- fread("../STAPSeq_vanja/db/normalized_counts/GSE116197_Enhancer_screens_Normalized_tagcounts_per_oligo.txt.gz")
DSCP <- DSCP[oligo_id=="DSCP"]
DSCP <- DSCP[, zfh1_enh_avg:= mean(c(zfh1_enh_Rep1, zfh1_enh_Rep2))]
lib <- rbind(lib, DSCP, fill= T)

#-------------------#
# Color data based on overlaps
#-------------------#
dat <- copy(lib)
dat[, basal:= log2(no_enh_Rep1+1)]
dat[, hk:= log2(ssp3_enh_Rep1+1)]
dat[, hk_norm:= hk-basal]
dat[, dev:= log2(zfh1_enh_avg+1)]
dat[, dev_norm:= dev-basal]
dat[, ratio:= hk-dev]

#-------------------#
# Color data based on overlaps
#-------------------#
# Franzis CP
dat[oligo_id=="chr3L_217307_217440_+", c("Cc", "name"):= .("red", "fl_8"), on= "oligo_id"]
# DSCP
dat[oligo_id=="DSCP", c("Cc", "name"):= .("limegreen", "DSCP")]
# Overlaps Rps12 promoter
RPS12_tss <- data.table(seqnames= "chr3L", 
                        start= 13016456, 
                        end= 13016456, 
                        strand= "+")
dat[RPS12_tss, c("Cc", "name"):= .("gold", "RpS12"), on= c("seqnames", "start<end", "end>start", "strand")]
#Check for the presence of restriction sites incompatible with STARR-Seq
dat[, enzyme_cut:= ifelse(grepl("ACCGGT|GTCGAC|CGCGCGCG", sequence), T, F)] # AgeI|SalI|MAUBI

#-------------------#
# Select new CPs
#-------------------#
optim_prim <- function(seq)
{
    .F <- rbindlist(lapply(0:4, function(x) c(strand= "F", seq= .c <- substr(seq, 1, 20+x), vl_oligo_Tm(.c))))
    .R <- rbindlist(lapply(0:4, function(x) c(strand= "R", seq= .c <- substr(vl_revComp(seq), 1, 20+x), vl_oligo_Tm(.c))))
    primers <- rbind(.F, .R)
    primers <- primers[between(`GC%`, 35, 65, incbounds = T) & between(Tm, 55, 58, incbounds = T)]
    return(primers[, .SD[1], strand])
}
cols <- c("name", "Cc", "F_primer", "R_primer")
dat[between(hk, 6.5, 8) & basal==0 & dev<1 & !enzyme_cut, (cols):= {
    prim <- optim_prim(sequence)
    if(nrow(prim)==2)
        .("vl_hk_low_", "blue", prim[1, seq], prim[2, seq]) else
            .(as.character(NA), as.character(NA), as.character(NA), as.character(NA))
    }, sequence]
dat[between(basal, 5, 5.5) & between(ratio, -7.1, -6.75) & dev>10.75 & !enzyme_cut, (cols):= {
    prim <- optim_prim(sequence)
    if(nrow(prim)==2)
        .("vl_dev_high_", "purple", prim[1, seq], prim[2, seq]) else
            .(as.character(NA), as.character(NA), as.character(NA), as.character(NA))
}, sequence]
dat[between(basal, 5.2, 5.5) & between(ratio, 2, 4) & !enzyme_cut, (cols):= {
    prim <- optim_prim(sequence)
    if(nrow(prim)==2)
        .("vl_hk_high_", "black", prim[1, seq], prim[2, seq]) else
            .(as.character(NA), as.character(NA), as.character(NA), as.character(NA))
}, sequence]
dat[grepl("vl_", name), name:= paste0(name, seq(.N)), name]

#-------------------#
# PLOT
#-------------------#
pdf("pdf/STARRSeq_design/RpS12_promoter_barplot.pdf", width= 10)
mat <- dat[Cc=="gold", .(oligo_id, basal, dev, hk)]
mat <- t(as.matrix(mat, 1))
barplot(mat, 
        beside= T, 
        las= 1,
        ylab= "log2FC")
legend("topright", 
       fill= gray.colors(nrow(mat)), 
       legend= c("basal", "dev", "hk"),
       bty= "n")
dev.off()

pdf("pdf/STARRSeq_design/CPs_hk_dev_responsiveness.pdf", width = 6, height = 6)
par(pty= "s", 
    las= 1, 
    mfrow= c(2,2), cex= 0.5)
# Plot hk
smoothScatter(dat$basal,
              dat$hk,
              colramp = colorRampPalette(c("white", "grey30", "grey20", "grey10")),
              xlab= "Basal activity",
              ylab= "hk")
points(dat[!is.na(Cc), basal],
       dat[!is.na(Cc), hk], 
       col= dat[!is.na(Cc), Cc],
       pch= 19)
text(dat[!is.na(Cc), basal],
     dat[!is.na(Cc), hk], 
     labels = dat[!is.na(Cc), name], 
     pos= 4,
     col= dat[!is.na(Cc), Cc])
# legend("bottomright", 
#        bty= "n",
#        pch= 19, 
#        col= c("red", "green", "gold", "blue"),
#        legend = c("Franzi", "vl", "RpS12", "AgeI/SalI"))
# Plot dev
smoothScatter(dat$basal,
              dat$dev,
              colramp = colorRampPalette(c("white", "grey30", "grey20", "grey10")),
              xlab= "Basal activity",
              ylab= "dev")
points(dat[!is.na(Cc), basal],
       dat[!is.na(Cc), dev], 
       col= dat[!is.na(Cc), Cc],
       pch= 19)
text(dat[!is.na(Cc), basal],
     dat[!is.na(Cc), dev], 
     labels = dat[!is.na(Cc), name], 
     pos= 4,
     col= dat[!is.na(Cc), Cc])
# Plot ratio
smoothScatter(dat$basal,
              dat$ratio,
              colramp = colorRampPalette(c("white", "grey30", "grey20", "grey10")),
              xlab= "Basal activity",
              ylab= "hk/dev")
points(dat[!is.na(Cc), basal],
       dat[!is.na(Cc), ratio], 
       col= dat[!is.na(Cc), Cc],
       pch= 19)
text(dat[!is.na(Cc), basal],
     dat[!is.na(Cc), ratio], 
     labels = dat[!is.na(Cc), name], 
     pos= 4,
     col= dat[!is.na(Cc), Cc])
# Plot hk vs dev
smoothScatter(dat$hk_norm,
              dat$dev_norm,
              colramp = colorRampPalette(c("white", "grey30", "grey20", "grey10")),
              xlab= "hk inducibility",
              ylab= "dev inducibility")
points(dat[!is.na(Cc), hk_norm],
       dat[!is.na(Cc), dev_norm], 
       col= dat[!is.na(Cc), Cc],
       pch= 19)
text(dat[!is.na(Cc), hk_norm],
     dat[!is.na(Cc), dev_norm], 
     labels = dat[!is.na(Cc), name], 
     pos= 4,
     col= dat[!is.na(Cc), Cc])
dev.off()

#------------------------------#
# SAVE selected primers
#------------------------------#
sel <- dat[!is.na(F_primer), .(oligo_id, name, seqnames, start, end, strand, F_primer, R_primer, CP_sequence= sequence)]
fwrite(sel,
       "Rdata/selected_CPs_primers_final.txt",
       sep= "\t")

#------------------------------#
# Check COFs inducibility
#------------------------------#
STAP <- get(load("/groups/stark/haberle/projects/Drosophila_STAPseq/results/factor.recruitment/CoF_screen_oligo_pool_selection_from_all_CoFs_for_first_submission_20171220/correlation.and.clustering/oligos.clustering/merged.reps/Oligos.clustering.based.on.Z-score.kmeans.5.of.total.tagcount.at.sig.activated.oligos.above.5.total.tags.and.3.tags.at.dominant.tss.in.at.least.two.reps.TSS.assigned.clusters.RData"))
STAP <- as.data.table(STAP)
cols <- names(STAP)[3:(ncol(STAP)-3)]
STAP[, paste0(cols, "_norm") := lapply(.SD, function(x) scale(log2(x+1)-log2(GFP+1))), .SDcols= cols]
all_CPs <- c(sel$oligo_id, 
             "chr3L_13016389_13016522_+", 
             "chr3L_13016403_13016536_+", 
             "chr3L_13016429_13016562_+", 
             "chr3L_217307_217440_+",
             "DSCP")
mat <- as.matrix(STAP[all_CPs, .SD, on= "bowtie_index", .SDcols= patterns("_norm")])
rownames(mat) <- paste0(all_CPs, "/", dat[all_CPs, name, on= "oligo_id"], " -> cl", STAP[sel$oligo_id, TSS_cluster, on= "bowtie_index"])

pdf("pdf/STARRSeq_design/heatmp_COF_induction_active_CPs.pdf", width = 12)
vl_heatmap(mat, 
           breaks = c(-2,0,8), 
           display_numbers = T)
dev.off()