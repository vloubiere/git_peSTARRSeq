
pdf("pdf/draft/screenshot_hk_vs_dev_STARRSeq.pdf", width = 6.5, height = 2.5)
par(mar= c(5,15,2,2))
regions <- data.table(seqnames= "chr3R", 
                      start= 5334856, 
                      end= 5344758)
regions[, {
  vl_screenshot(data.table(seqnames, start, end), 
                tracks = c("../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
                           "../gw_STARRSeq_bernardo/db/bw/RpS12_200bp_gw.UMI_cut_merged.bw"), 
                col= c("#74C37A", "tomato"), 
                genome = "dm3", 
                max = c(270, 170))
}, (regions)]

dev.off()

