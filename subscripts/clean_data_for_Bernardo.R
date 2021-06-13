feat <- readRDS("Rdata/uniq_library_final.rds")
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")

dat[feat, c("group_L", "detail_L", "coor_L"):= .(i.group, i.detail, paste0(i.seqnames, ":", i.start, "-", i.end)), on= "L==ID"]
dat <- dat[!is.na(group_L)]
dat[feat, c("group_R", "detail_R", "coor_R"):= .(i.group, i.detail, paste0(i.seqnames, ":", i.start, "-", i.end)), on= "R==ID"]
dat <- dat[!is.na(group_R)]

final <- dat[, .(L, R, coor_L, coor_R, group_L, group_R, detail_L, detail_R, 
                 act_L= median_L, act_R= median_R, act_pair= log2FoldChange, additive_score= add, residuals= log2FoldChange-add, lib)]
final[, {
  saveRDS(.SD, paste0("/groups/stark/vloubiere/projects/pe_STARRSeq_2/for_Bernado/", lib, "_data.rds"))
}, lib]


# Dear Bernardino, 
# Here are the links for my 3 different libraries:
# - Small (10,000 pairs): "/groups/stark/vloubiere/projects/pe_STARRSeq_2for_Bernado/libvl013_data.rds"
# - Medium (160,000 pairs): "/groups/stark/vloubiere/projects/pe_STARRSeq_2for_Bernado/libvl014_data.rds"
# - Big (1,000,000 pairs): "/groups/stark/vloubiere/projects/pe_STARRSeq_2for_Bernado/libvl002_data.rds"
# 
# They all contain the same columns:
# - L: Left enhancer ID
# - R: Right enhancer ID
# - coor_L: Left enhancer genomic coordinates 
# - coor_R: Right enhancer genomic coordinates 
# - group_L: Left enhancer group ("dev", "hk", "inducible" (hs/ecd), "OSC" (OSC-specific), "shared" (hk+dev))
# - group_R: Right enhancer group ("dev", "hk", "inducible" (hs/ecd), "OSC" (OSC-specific), "shared" (hk+dev))
# - detail_L: Left enhancer detail ("inactive", "weak", "medium", "strong", "ecdysone", "heatshock")
# - detail_R: Right enhancer detail ("inactive", "weak", "medium", "strong", "ecdysone", "heatshock")
# - act_L: Left enhancer individual activity (log2)
# - act_R: Right enhancer individual activity (log2)
# - act_pair: Pair activity (log2)
# - additive_score: Pair expected additive activity (log2)
# - residuals: Cooperativity residuals (= act_pair-additive_score)
# 
# For now, I could not find any sequence feature that would be specific of the most super-additive pairs, and this is the main point I would like you to double check. The only thing I found so far in terms of sequence features, is that the DREF motif is negatively associated with super-additivity (consistent with the fact that I use the DSCP developmental promoter)
#
# Finally, a few technical aspects:
# 1/ You could start with the smaller library ("vllib013"), made of ~10,000 pairs in total. However, you will notice that you only get ~6000 lines which is a consequence of 1/ few missing pairs and 2/ the fact that control containing-pairs have been removed during late processing steps.
# 2/ When doing my analyses, I only consider the pairs for which both enhancers are individually active (act_L>0.5 & act_R>0.5).
# 3/ Also, be careful about the pairs containing >=1 hk enhancer (group_L=="hk" | group_R=="hk"), which behave differently.