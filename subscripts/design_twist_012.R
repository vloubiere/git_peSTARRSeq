setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(BSgenome.Dmelanogaster.UCSC.dm3)
peSTARR <- fread("old_versions/20210613_ce87ade_backup/db/DE_analysis/libvl002_DE.txt")
dat <- as.data.table(readRDS("old_versions/20210613_ce87ade_backup/Rdata/master_results_peSTARRSeq.rds"))
lib <- as.data.table(readRDS("old_versions/20210613_ce87ade_backup/Rdata/uniq_library_final.rds"))
lib_seq <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
lib[lib_seq, enh_seq:= i.enh_sequence, on= "ID==ID_vl"]

#-----------------------------#
# Negative controls
#-----------------------------#
ctls <- lib[group=="control" & (vl)]
# Identify controls which give similar results regardless of position
ctls[, PCC:= {
  .c <- merge(peSTARR[.BY, , on= "L==ID"],
              peSTARR[.BY, , on= "R==ID"],
              by.x= "R",
              by.y= "L")
  if(nrow(.c)>10)
    cor.test(.c$log2FoldChange.x, .c$log2FoldChange.y)$estimate else
      as.numeric(NA)
}, ID]
# Identify controls which have truely low/null basal activity
peSTARR[grepl("control", L), out_L:= {
  .c <- boxplot(log2FoldChange, plot= F)
  !between(log2FoldChange, .c$stats[1,1], .c$stats[5,1])
}, R]
peSTARR[grepl("control", R), out_R:= {
  .c <- boxplot(log2FoldChange, plot= F)
  !between(log2FoldChange, .c$stats[1,1], .c$stats[5,1])
}, L]
ctls[, c("outliers_freq_L", "outliers_freq_R"):= {
  .L <- unique(peSTARR[.BY, , on= "L==ID"][, .(out_L, R)])
  .R <- unique(peSTARR[.BY, , on= "R==ID"][, .(out_R, L)])
  .(length(which(.L$out_L))/nrow(.L),
    length(which(.R$out_R))/nrow(.R))
}, ID]
ctl_sel <- ctls[PCC>0.8 & outliers_freq_L<0.025 & outliers_freq_R<0.025, ID]
ctl_sel <- lib[ID %in% ctl_sel]
ctl_sel <- ctl_sel[, .(seqnames, start, end, strand, group, detail, enh_seq)]

#-----------------------------#
# Dev enhancers
#-----------------------------#
L_act <- unique(dat[lib=="libvl014" & median_L>1 & grepl("dev", L), L])
R_act <- unique(dat[lib=="libvl014" & median_R>1 & grepl("dev", R), R])
dev_sel <- L_act[L_act %in% R_act & grepl("_B_", L_act)]
dev_sel <- lib[ID %in% dev_sel]
dev_sel <- dev_sel[, .(seqnames, start, end, strand, group, detail, enh_seq)]

#-----------------------------#
# hk enhancers
#-----------------------------#
hk_sel <- lib[(vl) & group=="hk" & !grepl("inactive", ID), ID]
set.seed(1)
hk_sel <- lib[ID %in% sample(hk_sel, nrow(dev_sel))]
hk_sel <- hk_sel[, .(seqnames, start, end, strand, group, detail, enh_seq)]

#-----------------------------#
# SUHW peaks
#-----------------------------#
suhw_sel <- readRDS("Rdata/SUHW_peaks.rds")
suhw_sel <- suhw_sel[, .(seqnames, 
                         start= max_coor-124, 
                         end= max_coor+124, 
                         strand= "*",
                         group= "SUHW_peak",
                         detail= "SUHW_peak")]
suhw_sel[, enh_seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(suhw_sel)))]

#-----------------------------#
# DHS+ STARR- peaks
#-----------------------------#
DHS_sel <- as.data.table(readRDS("Rdata/DHS+_STARR-_sequences.rds"))
DHS_sel <- DHS_sel[, .(seqnames,
                       start= round(rowMeans(.SD))-124, 
                       end= round(rowMeans(.SD))+124, 
                       strand= "*",
                       group= "DHS_peak",
                       detail= "no_STARR"), .SDcols= c("start", "end")]
DHS_sel[, enh_seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(DHS_sel)))]

#-----------------------------#
# Repressors
#-----------------------------#
rep_sel <- fread("Rdata/Lorena_top_repressors_2010_jung_DLM3+521_twist_oligos_VL.txt")
rep_sel <- rep_sel[, .(seqnames, 
                       start= round(rowMeans(.SD))-124, 
                       end= round(rowMeans(.SD))+124, 
                       strand= "*",
                       group= "Repressor",
                       detail= "Top_dip_LH"), .SDcols= c("start", "end")]
rep_sel[, enh_seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(rep_sel)))]

#-----------------------------#
# CPs
#-----------------------------#
CP_sel <- readRDS("Rdata/CPs_selection_twist.rds")
CP_sel[strand=="+", end:= end+44] # Discussed with Vanja (orginally +66 downstream)
CP_sel[strand=="+", start:= end-248]
CP_sel[strand=="-", start:= start-44] # Discussed with Vanja (orginally +66 downstream)
CP_sel[strand=="-", end:= start+248]
CP_sel <- CP_sel[, .(seqnames, 
                     start= start, 
                     end= end, 
                     strand,
                     group= "CP",
                     detail)]
CP_sel[, enh_seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(CP_sel[, .(seqnames, start, end)])))]

#-----------------------------#
# Final_object
#-----------------------------#
final <- rbind(ctl_sel, 
               dev_sel, 
               hk_sel, 
               CP_sel, 
               suhw_sel, 
               DHS_sel, 
               rep_sel)
linkers <- unique(lib_seq[, .(linker_ID, fw_linker, rev_linker)])
setkeyv(linkers, "linker_ID")
final[group %in% c("control", "dev", "hk", "shared"), c("fw_linker", "rev_linker", "linker_ID"):= .(linkers["A"]$fw_linker, 
                                                                                                    linkers["A"]$rev_linker,
                                                                                                    "A")]
final[group %in% c("DHS_peak", "Repressor"), c("fw_linker", "rev_linker", "linker_ID"):= .(linkers["B"]$fw_linker, 
                                                                                           linkers["B"]$rev_linker,
                                                                                           "B")]
final[group %in% c("CP", "SUHW_peak"), c("fw_linker", "rev_linker", "linker_ID"):= .(linkers["C"]$fw_linker, 
                                                                                     linkers["C"]$rev_linker,
                                                                                     "C")]
final[, oligo_full_sequence:= paste0(fw_linker, enh_seq, rev_linker)]
setorderv(final, c("group", "detail", "seqnames", "start"))
final[, ID:= paste0(group,  
                    ifelse(detail %in% c("inactive", "weak", "medium", "strong"), paste0("_", detail), ""))]
final[, ID:= paste0(ID, "_", linker_ID, "_"), linker_ID]
final[, ID:= paste0(ID, sprintf("%05d", .I))]

saveRDS(final, "Rdata/vl_library_twist12_210610.rds")






