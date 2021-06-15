setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
avail_dir <- "/groups/stark/vloubiere/projects/available_data_dm3/"
require(vlfunctions)
require(BSgenome.Dmelanogaster.UCSC.dm3)

#---------------------------------------#
# ATAC-Seq peak calling
#---------------------------------------#
if(!file.exists("db/peaks/ATAC_peaks.txt"))
{
  peaks <- vl_peakCalling(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep1_uniq.bed",
                                       "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep2_uniq.bed",
                                       "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep3_uniq.bed",
                                       "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep4_uniq.bed"))
  fwrite(peaks[OR>2], 
         "db/peaks/ATAC_peaks.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}

#---------------------------------------#
# DHS peak calling
#---------------------------------------#
if(!file.exists("db/peaks/DHS_peaks.txt"))
{
  peaks <- vl_peakCalling(ChIP_bed = "/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/mapped//S2_DHS.mapped.reads.bed",
                          Input_bed = "/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/mapped//S2_DHS_input.mapped.reads.bed")
  fwrite(peaks[OR>3], 
         "db/peaks/DHS_peaks.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}

#---------------------------------------#
# STARR-Seq peak calling
#---------------------------------------#
if(!file.exists("db/peaks/STARR_DSCP_200_peaks.txt"))
{
  peaks <- vl_peakCalling(ChIP_bed = c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed",
                                       "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed"),
                          Input_bed =  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/input_DSCP_200bp_cut.bed")
  fwrite(peaks[OR>3], 
         "db/peaks/STARR_DSCP_200_peaks.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}
if(!file.exists("db/peaks/STARR_DSCP_600_peaks.txt"))
{
  peaks <- vl_peakCalling(ChIP_bed = c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data//DSCP_600bp_gw_Rep1_cut.bed" ,
                                       "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data//DSCP_600bp_gw_Rep2_cut.bed" ),
                          Input_bed =  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data//input_DSCP_600bp_cut.bed")
  fwrite(peaks[OR>4], 
         "db/peaks/STARR_DSCP_600_peaks.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}

#---------------------------------------#
# H3K4me3 peak calling
#---------------------------------------#
if(!file.exists("db/peaks/STARR_DSCP_200_peaks.txt"))
{
  peaks <- vl_peakCalling(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep1_uniq.bed",
                                       "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep2_uniq.bed"),
                          Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep1_uniq.bed",
                                        "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep2_uniq.bed"))
  fwrite(peaks, 
         "db/peaks/H3K4me3_peaks.txt", 
         col.names = T, 
         row.names = F, 
         sep= "\t", 
         quote= F, 
         na= NA)
}