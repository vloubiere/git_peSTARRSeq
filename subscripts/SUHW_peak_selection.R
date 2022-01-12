setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(BSgenome.Dmelanogaster.UCSC.dm3)

#--------------------#
# Peak calling
#--------------------#
if(!file.exists("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_SuHw_rep1_uniq.bed"))
{
  SUHW <- vl_peakCalling("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41354_SuHw_rep1_uniq.bed")
  
  fwrite(SUHW[OR>4], 
         "db/peaks/SUHW_peaks.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}else
{
  SUHW <- fread("db/peaks/SUHW_peaks.txt")
}

#--------------------#
# Measure STARR-Seq enrichment
#--------------------#
SUHW$STARR200 <- vl_computeEnrichment(ChIP_bed = c("../gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed",
                                                   "../gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed"),
                                      Input_bed = "../gw_STARRSeq_bernardo/db/raw_data/input_DSCP_200bp_cut.bed",
                                      peaks =  SUHW,
                                      ext_peaks = 50000)$log2_enr

SUHW$STARR600 <- vl_computeEnrichment(ChIP_bed = c("../gw_STARRSeq_bernardo/db/raw_data//DSCP_600bp_gw_Rep1_cut.bed" ,
                                                   "../gw_STARRSeq_bernardo/db/raw_data//DSCP_600bp_gw_Rep2_cut.bed" ),
                                      Input_bed = "../gw_STARRSeq_bernardo/db/raw_data//input_DSCP_600bp_cut.bed",
                                      peaks =  SUHW,
                                      ext_peaks = 50000)$log2_enr

# plot(SUHW$STARR200, SUHW$OR)
# sel <- identify(SUHW$STARR200, SUHW$OR)
sel <- c(36, 102, 488, 608, 1044, 1131, 1266, 1336, 1345, 1600, 1670, 1673, 1796, 2520, 2625, 3065, 3101, 3146, 3242, 3251,
         3328, 3346, 3575, 3587, 3755, 4098, 4494, 4546, 4563, 4577, 4642, 4687, 4864, 4899, 4926, 4948, 4955, 4961, 5026, 
         5028, 5033, 5093, 5180, 5193, 5281, 5293, 5318, 5477, 5527, 5557, 5582, 404, 898, 1149, 1451, 1965, 2030, 2449, 
         2963, 3138, 3308, 3310, 3486, 3838, 3992, 4435, 4647, 4662, 4831, 5150, 5213)
plot(SUHW$STARR200, 
     SUHW$OR, 
     col= ifelse(seq(nrow(SUHW)) %in% sel, "red", adjustcolor("lightgrey", 0.3)), 
     pch= 19)

saveRDS(SUHW[sel], "Rdata/SUHW_peaks.rds")

