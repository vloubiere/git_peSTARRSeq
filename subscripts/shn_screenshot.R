setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

region <- data.table(seqnames= "chr2R",
                     start= 7046591,
                     end= 7104899)
tracks <- c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw",
            "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE119708_ATAC_rep2_uniq.bw",
            "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41440_H3K27ac_rep1_uniq.bw")

dir.create("pdf/screnshots", showWarnings = F)

pdf("pdf/screnshot_shn_STARRSeq.pdf", height = 4)
vl_screenshot(region, 
              tracks, 
              n_genes = 2, 
              gband = 0.3, 
              names = c("STARR-Seq", "ATAC-Seq", "H3K27Ac"))
dev.off()