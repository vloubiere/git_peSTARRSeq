setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir.create("pdf/average_tracks/", 
           showWarnings = F)

#-----------------------#
# IMPORT
#-----------------------#
dat <- readRDS("Rdata/final_300bp_enhancer_features.rds")

#-----------------------#
# SUHW peaks
#-----------------------#
vl_average_bw_track(bed = dat[group %in% c("Silencer", "SUHW_peak"), .(seqnames, start, end)], 
                    set_IDs = dat[group %in% c("Silencer", "SUHW_peak"), group], 
                    tracks = "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw")

#-----------------------#
# DHS peaks
#-----------------------#
tracks <- data.table(file= c("/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/tracks/S2_DHS.bw",
                             "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                             "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw"))

pdf("pdf/average_tracks/DHS_peaks_vs_dev_enhancers.pdf")
par(mfrow= c(3, 1),
    pty= "s")
tracks[, {
  vl_average_bw_track(bed = dat[group %in% c("DHS_peak", "dev"), .(seqnames, start, end)], 
                      set_IDs = dat[group %in% c("DHS_peak", "dev"), group], 
                      extend = c(-2000, 2000),
                      tracks = file)
  print("")
}, file]
dev.off()

