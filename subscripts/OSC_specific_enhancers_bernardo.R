setwd("~/Dropbox (VBC)/R_db/dm3/0002_library_design/")

# load OSC enhancer summits
load("../Bernardo_fanny/Rdata/All_Drosophila_STARRseq_peaks_and_negative.RData")
OSC_enhancers <- Drosophila_enhancers[Drosophila_enhancers$Group %in% "OSC"]

# remove the ones overlapping S2 enhancers
OSC_specific_enhancers <- OSC_enhancers[-queryHits(findOverlaps(resize(OSC_enhancers,501,"center"),
                                                                resize(Drosophila_enhancers[Drosophila_enhancers$Group %in% c("S2_Science2013", "dCP", "hkCP")],501,"center")))]

# resize to 249bp
OSC_specific_enhancers <- resize(OSC_specific_enhancers, 249, "center")
strand(OSC_specific_enhancers) <- "+"

# load calculated enrichments in S2 STARR-seq (dev and hk, 600bp and 200bp) - to make sure I select OSC enhancers with no activity in S2 cells
Test_S2_enrichments <- readRDS("../Bernardo_fanny/Rdata/OSC_and_neg_tested_S2_enrichments.rds")

OSC_specific_enhancers <- OSC_specific_enhancers[Test_S2_enrichments$OSC$hk_600bp_log2_enr <=0 &
                                                   Test_S2_enrichments$OSC$dev_600bp_log2_enr <=0 &
                                                   Test_S2_enrichments$OSC$hk_200bp_log2_enr <=0 &
                                                   Test_S2_enrichments$OSC$dev_200bp_log2_enr <=0] # because order is the same, I don't need to match
OSC_specific_enhancers <- OSC_specific_enhancers[order(OSC_specific_enhancers$enrichment, decreasing = T)]
save(OSC_specific_enhancers, file= "../Bernardo_fanny/Rdata/OSC_specific_enhancers.Rdata")
