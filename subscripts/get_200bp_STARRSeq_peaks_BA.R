load_peak_function <- function(x, width){
  path<- "/groups/stark/almeida/data/STARRseq/dm3/20190401_200bp_gw_STARRseq/data/"
  # load peaks file and convert to granges
  gr <- makeGRangesFromDataFrame(read.delim(paste0(path, x), header = F),
                                 seqnames.field = "V1", start.field = "V2", end.field = "V2", keep.extra.columns = T)
  # resize peaks
  gr <- resize(gr, width, "center")
  # Treat peak information
  names(mcols(gr))[5:7] <- c("Enrch.", "Corr_enrch", "p_value")
  gr
}



peaks_list <- list(dCP_200bp = load_peak_function("DSCP_200bp_gw.UMI_cut_merged.peaks.txt", 201),
                   hkCP_200bp = load_peak_function("RpS12_200bp_gw.UMI_cut_merged.peaks.txt", 201))
