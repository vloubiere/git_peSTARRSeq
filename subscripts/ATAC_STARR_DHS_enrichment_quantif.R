setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#---------------------------#
# Quantif signal
#---------------------------#
peaks <- data.table(file= c("db/peaks/ATAC_peaks.txt",
                            "db/peaks/DHS_peaks.txt",
                            "db/peaks/H3K4me3_peaks.txt",
                            "db/peaks/STARR_DSCP_200_peaks.txt",
                            "db/peaks/STARR_DSCP_600_peaks.txt",
                            "db/peaks/SUHW_peaks.txt"))
peaks <- peaks[, fread(file, colClasses = c("character", rep("numeric", 5))), file]
peaks <- resize(GRanges(peaks$seqnames, IRanges(peaks$max_coor, peaks$max_coor)), 100, "center")
peaks <- vl_collapse_DT_ranges(as.data.table(peaks))
peaks <- unique(as.data.table(resize(GRanges(peaks), width = 500, "center")))

peaks$H3K4me1 <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep1_uniq.bed",
                                                   "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K4me1_rep2_uniq.bed"),
                                      Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                    "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"),
                                      peaks =  peaks,
                                      ext_peaks = 50000)$log2_enr

peaks$H3K4me3 <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep1_uniq.bed",
                                                   "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_H3K4me3_rep2_uniq.bed"),
                                      Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep1_uniq.bed",
                                                    "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE81795_input_rep2_uniq.bed"),
                                      peaks =  peaks,
                                      ext_peaks = 50000)$log2_enr

peaks$H3K27ac <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep1_uniq.bed",
                                                   "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27ac_rep2_uniq.bed"),
                                      Input_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                                    "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep2_uniq.bed"),
                                      peaks =  peaks,
                                      ext_peaks = 50000)$log2_enr

peaks$H3K27me3 <- vl_computeEnrichment(ChIP_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_H3K27me3_rep1_uniq.bed",
                                       Input_bed = "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE41440_input_rep1_uniq.bed",
                                       peaks =  peaks,
                                       ext_peaks = 50000)$log2_enr

peaks$ATAC <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep1_uniq.bed",
                                                "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep2_uniq.bed",
                                                "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep3_uniq.bed",
                                                "/groups/stark/vloubiere/projects/available_data_dm3/db/bed/GSE119708_ATAC_rep4_uniq.bed"),
                                   peaks =  peaks,
                                   ext_peaks = 50000)$log2_enr

peaks$DHS <- vl_computeEnrichment(ChIP_bed = "/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/mapped//S2_DHS.mapped.reads.bed",
                                  Input_bed = "/groups/stark/haberle/data/public/dm3/S2_cells/DHS/Stark_S2_cells/mapped//S2_DHS_input.mapped.reads.bed",
                                  peaks =  peaks,
                                  ext_peaks = 50000)$log2_enr

peaks$STARR200 <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed",
                                                    "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed"),
                                       Input_bed =  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/input_DSCP_200bp_cut.bed",
                                       peaks =  peaks,
                                       ext_peaks = 50000)$log2_enr

peaks$STARR600 <- vl_computeEnrichment(ChIP_bed = c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data//DSCP_600bp_gw_Rep1_cut.bed" ,
                                                    "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data//DSCP_600bp_gw_Rep2_cut.bed" ),
                                       Input_bed =  "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data//input_DSCP_600bp_cut.bed",
                                       peaks =  peaks,
                                       ext_peaks = 50000)$log2_enr

fwrite(peaks, 
       file = "db/library_design/twist012/STARR_ATAC_DHS_CHIP_quantif.txt", 
       col.names = T, 
       row.names = F, 
       sep= "\t", 
       quote= F,
       na= NA)



