dat <- fread("old_versions/original_folder_twist12_design_backup/Rdata/STARR_ATAC_DHS_CHIP_quantif.txt")
dat[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
K4me3 <- fread("old_versions/original_folder_twist12_design_backup/db/peaks/H3K4me3_peaks.txt")
K4me3[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
STARR200 <- fread("old_versions/original_folder_twist12_design_backup/db/peaks/STARR_DSCP_200_peaks.txt")
STARR200[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
STARR600 <- fread("old_versions/original_folder_twist12_design_backup/db/peaks/STARR_DSCP_600_peaks.txt")
STARR600[, center:= rowMeans(.SD), .SDcols= c("start", "end")]

dat$K4me3_dist <- K4me3[dat, min(abs(center-i.center)), .EACHI, on= "seqnames"]$V1
dat$STARR200_dist <- STARR200[dat, min(abs(center-i.center)), .EACHI, on= "seqnames"]$V1
dat$STARR600_dist <- STARR600[dat, min(abs(center-i.center)), .EACHI, on= "seqnames"]$V1

sel <- dat[(ATAC>2.5 | DHS>2.5)
           & H3K27ac>1
           & H3K4me1>1
           & K4me3_dist>2500
           & STARR200_dist>2500
           & STARR600_dist>2500
           & seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")]

sel[, ext:= as.character(resize(GRanges(sel), 5000, "center"))]
sel[, peaks:= as.character(resize(GRanges(sel), 250, "center"))]
sel[, idx:= rep(seq(1000), each= 5)[1:.N]]

bw <- c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
        "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw",
        "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE119708_ATAC_rep1_uniq.bw",
        "/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/tracks/S2_DHS.bw",
        # "/groups/stark/vlasova/data/STARRseq/dm3/Lorena_repression_merged/samples/libLOH009_enh10-2_DSCP_gw200_rep_1.UMI_ps.bw",
        # "/groups/stark/vlasova/data/STARRseq/dm3/Lorena_repression_merged/samples/libLOH009_enh10-2_DSCP_gw200_rep_1.UMI_ns.bw",
        # "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSM2055132_S2_NHS_PRO_seq_norm_plus.bw",
        # "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSM2055132_S2_NHS_PRO_seq_norm_minus.bw",
        # "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw",
        # "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41440_H3K27me3_rep1_uniq.bw"
        "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41440_H3K4me1_rep1_uniq.bw",
        "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw",
        "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41440_H3K27ac_rep1_uniq.bw")

pdf("pdf/screenshots_selection_DHS+_STARR-.pdf", width = 10, height = 7)
par(mar= c(1,15,1,1))
sel[, {
  vl_screenshot(GRanges(ext),
                bw,
                highlight_regions = GRanges(peaks),
                max = c(45, 100, 8, 8, 13, 35, 30))
  print("")
}, idx]
dev.off()

saveRDS(sel, "Rdata/DHS+_STARR-_sequences.rds") 

