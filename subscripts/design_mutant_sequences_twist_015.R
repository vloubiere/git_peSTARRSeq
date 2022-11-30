setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Clean TWIST library
if(!file.exists("db/library_design/twist015/candidates_clean_twist015_design.txt"))
{
  # import TWIST data
  TWIST <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
  lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
  TWIST[lib, vlID:= i.ID_vl, on= "ID==BA_ID"]
  dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")
  
  # Dev cutoff Left and Right enhancers (avoid saturation)
  cutoffL <- mean(TWIST[unique(dat[indL<4, .(L, indL)])[order(-indL), L][1:10], dev_log2FoldChange, on= "vlID"])
  cutoffR <- mean(TWIST[unique(dat[indR<4, .(R, indR)])[order(-indR), R][1:10], dev_log2FoldChange, on= "vlID"])
  TWIST <- TWIST[dev_log2FoldChange < min(c(cutoffL, cutoffR))]
  TWIST[, enh_sequence:= vl_getSequence(TWIST, "dm3")]
  
  # Import DHS Data
  DHS <- fread("/groups/stark/almeida/Projects/Sequence_rules_epigenomic_layers/results/20220125_Enh_act_and_chromatin/DHS_peaks_annotated_enh_act.txt")
  DHS[, summit:= start+summit]
  DHS[, start:= summit-124]
  DHS[, end:= summit+124]
  DHS[, enh_sequence:= vl_getSequence(DHS, "dm3")]
  DHS <- DHS[Enh_short_overlap=="Inactive"]
  
  # Merge
  lib <- rbind(TWIST, DHS, fill= T)
  lib <- lib[, .(BAID= ID, 
                 seqnames, start, end, strand, 
                 dev_log2FoldChange, dev_padj, hk_log2FoldChange, hk_padj, 
                 enh_sequence)]
  lib[is.na(BAID), BAID:= paste0("DHS_", .I)]
  
  # Motif positions
  pos <- vl_motif_pos(lib$enh_sequence, 
                      sel= c("flyfactorsurvey__CG16778_SANGER_5_FBgn0003715", 
                             "cisbp__M2014",
                             "homer__AVYTATCGATAD_DREF"),
                      genome= "dm3", 
                      collapse_overlapping = T,
                      p.cutoff = 1e-4)
  pos <- as.data.table(pos)
  setnames(pos, c("motif_twist_pos", "motif_Trl_pos", "motif_Dref_pos"))
  counts <- as.data.table(lapply(pos, function(x) sapply(x, nrow)))
  setnames(counts, c("motif_twist", "motif_Trl", "motif_Dref"))
  lib <- cbind(lib, counts, pos)
  saveRDS(lib, 
          "db/library_design/twist015/candidates_clean_twist015_design.txt")
}else
  lib <- readRDS("db/library_design/twist015/candidates_clean_twist015_design.txt")

#----------------------------------#
# WT sequences
#----------------------------------#
WT <- rbindlist(list(control= lib[hk_padj>0.05 & dev_padj>0.05 & dev_log2FoldChange<1 & hk_log2FoldChange<1 
                                  & motif_twist+motif_Trl+motif_Dref==0
                                  & !grepl("^DHS", BAID)], 
                     DHS= lib[grepl("^DHS", BAID) & motif_twist+motif_Trl==0], 
                     noMotifAct= lib[motif_twist+motif_Trl+motif_Dref==0 & dev_log2FoldChange>1], 
                     twist= lib[between(motif_twist, 2, 5)], # Avoid weird sequences full of GA
                     Trl= lib[between(motif_Trl, 3, 6)], # Avoid weird sequences full of GA
                     shared= lib[between(motif_Dref, 2, 5) & dev_log2FoldChange>1 & hk_log2FoldChange>1]),
                idcol = "group")
WT[group=="twist", motif:= lapply(motif_twist_pos, function(x) x[, c((start-8:end+8)), .(start, end)]$V1)]
WT[group=="twist", motif_center:= lapply(motif_twist_pos, function(x) round(rowMeans(x[, start:end])))]
WT[group=="Trl", motif:= lapply(motif_Trl_pos, function(x) x[, c(start-8:end+8), .(start, end)]$V1)]
WT[group=="Trl", motif_center:= lapply(motif_Trl_pos, function(x) round(rowMeans(x[, start:end])))]
WT[group=="shared", motif:= lapply(motif_Dref_pos, function(x) x[, c(start-8:end+8), .(start, end)]$V1)]
WT[group=="shared", motif_center:= lapply(motif_Dref_pos, function(x) round(rowMeans(x[, start:end])))]
# Unique sequences
WT <- WT[, .SD[1], BAID]

#------------------------------#
# Paste x2 Trl/twist Motifs in control/DHS/noMotifAct sequences
#------------------------------#
addSynMot <- CJ(pos1= 20:110, pos2= 140:230)
set.seed(1)
addSynMot <- addSynMot[sample(nrow(addSynMot), 1000)]
addSynMot <- WT[group %in% c("control", "DHS", "noMotifAct"), as.list(addSynMot), .(BAID, group, enh_sequence)]
addSynMot[, seq1:= substr(enh_sequence, 1, pos1-1), pos1]
addSynMot[, seq2:= substr(enh_sequence, pos1+8, pos2-1), .(pos1, pos2)]
addSynMot[, seq3:= substr(enh_sequence, pos2+8, 249), pos2]
# Trl
set.seed(1)
mot1 <- sample(c("GAGAGAGA", "CTCTCTCT"), nrow(addSynMot), replace = T)
set.seed(2)
mot2 <- sample(c("GAGAGAGA", "CTCTCTCT"), nrow(addSynMot), replace = T)
addSynMot$addTrl <- paste0(addSynMot$seq1, mot1, addSynMot$seq2, mot2, addSynMot$seq3)
# twist
set.seed(1)
mot1 <- sample(c("CACATATG", "CATATGTG"), nrow(addSynMot), replace = T)
set.seed(2)
mot2 <- sample(c("CACATATG", "CATATGTG"), nrow(addSynMot), replace = T)
addSynMot$addTwist <- paste0(addSynMot$seq1, mot1, addSynMot$seq2, mot2, addSynMot$seq3)

#------------------------------#
# Paste x3 Trl/twist Motifs in control/DHS sequences
#------------------------------#
addSynMotInact <- CJ(pos1= 20:110, pos2= 30:220, pos3= 140:230)
addSynMotInact <- addSynMotInact[pos2>pos1+10 & pos3>pos2+10]
set.seed(1)
addSynMotInact <- addSynMotInact[sample(nrow(addSynMotInact), 1000)]
addSynMotInact <- WT[group %in% c("control", "DHS"), as.list(addSynMotInact), .(BAID, group, enh_sequence)]
addSynMotInact[, seq1:= substr(enh_sequence, 1, pos1-1), pos1]
addSynMotInact[, seq2:= substr(enh_sequence, pos1+8, pos2-1), .(pos1, pos2)]
addSynMotInact[, seq3:= substr(enh_sequence, pos2+8, pos3-1), .(pos2, pos3)]
addSynMotInact[, seq4:= substr(enh_sequence, pos3+8, 249), pos3]
# Trl
set.seed(1)
mot1 <- sample(c("GAGAGAGA", "CTCTCTCT"), nrow(addSynMotInact), replace = T)
set.seed(2)
mot2 <- sample(c("GAGAGAGA", "CTCTCTCT"), nrow(addSynMotInact), replace = T)
set.seed(3)
mot3 <- sample(c("GAGAGAGA", "CTCTCTCT"), nrow(addSynMotInact), replace = T)
addSynMotInact$addTrl <- paste0(addSynMotInact$seq1, mot1, addSynMotInact$seq2, mot2, addSynMotInact$seq3, mot3, addSynMotInact$seq4)
# twist
set.seed(1)
mot1 <- sample(c("CACATATG", "CATATGTG"), nrow(addSynMotInact), replace = T)
set.seed(2)
mot2 <- sample(c("CACATATG", "CATATGTG"), nrow(addSynMotInact), replace = T)
set.seed(3)
mot3 <- sample(c("CACATATG", "CATATGTG"), nrow(addSynMotInact), replace = T)
addSynMotInact$addTwist <- paste0(addSynMotInact$seq1, mot1, addSynMotInact$seq2, mot2, addSynMotInact$seq3, mot3, addSynMotInact$seq4)

#------------------------------#
# Paste Dref Motifs in active/synergistic sequences
#------------------------------#
addDrefMot <- CJ(pos1= 20:110, pos2= 140:230)
set.seed(1)
addDrefMot <- addDrefMot[sample(nrow(addDrefMot), 1000)]
addDrefMot <- WT[group %in% c("noMotifAct", "Trl", "twist"), as.list(addDrefMot), .(BAID, group, enh_sequence)]
addDrefMot[, seq1:= substr(enh_sequence, 1, pos1-1), pos1]
addDrefMot[, seq2:= substr(enh_sequence, pos1+8, pos2-1), .(pos1, pos2)]
addDrefMot[, seq3:= substr(enh_sequence, pos2+8, 249), pos2]
# Remove positions overlapping Trl/twist motifs
rm <- WT[group %in% c("Trl", "twist"), motif[[1]], BAID]
addDrefMot[, keep:= TRUE]
addDrefMot[rm, keep:= FALSE, on= c("BAID", "pos1==V1")]
addDrefMot[rm, keep:= FALSE, on= c("BAID", "pos2==V1")]
addDrefMot <- addDrefMot[(keep)]
# Paste motif
addDrefMot$addDref <- paste0(addDrefMot$seq1, "TATCGATA", addDrefMot$seq2, "TATCGATA", addDrefMot$seq3)

#------------------------------#
# Mutate Trl/twist in synergistic enhancers
#------------------------------#
mutSynMot <- WT[group %in% c("Trl", "twist"), .(pos= motif_center[[1]]), .(BAID, group, enh_sequence)]
mutSynMot <- mutSynMot[, {
  exp <- if(group=="twist")
    rep(list(c("TAAGA", "TTGACC", "AAGTA", "CCTTAAG", "GCTTAA", "TCTTA", "GGTCAA", "TACTT", "CTTAAGG", "TTAAGC")), 
        length(pos)) else if(group=="Trl")
                 rep(list(c("TAAGA", "TTGACC", "AAGTA", "CCTTAAG", "GCTTAA")), length(pos))
  exp <- c(quote(CJ), exp)
  cmb <- eval(as.call(exp))
  cmb[, idx:= .I]
  cmb <- cmb[, .(pos, mot= unlist(.SD)), .(idx)]
  res <- cmb[, {
    .c <- enh_sequence
    for(i in seq(.N))
    {
      motWidth <- nchar(mot[i])
      is <- pos[i]-round(motWidth/2)
      ie <- is+motWidth-1
      substr(.c, is, ie) <- mot[i]
    }
    .c
  }, idx]$V1
  .(mutSyn= unique(res))
}, .(BAID, group)]

#------------------------------#
# Mutate Dref in shared enhancers
#------------------------------#
mutDrefMot <- WT[group=="shared", .(pos= motif_center[[1]]), .(BAID, group, enh_sequence)]
mutDrefMot <- mutDrefMot[, {
  exp <- rep(list(c("TAAGA", "TTGACC", "AAGTA", "CCTTAAG", "GCTTAA", "TCTTA", "GGTCAA", "TACTT", "CTTAAGG", "TTAAGC")), 
             length(pos))
  exp <- c(quote(CJ), exp)
  cmb <- eval(as.call(exp))
  cmb[, idx:= .I]
  cmb <- cmb[, .(pos, mot= unlist(.SD)), .(idx)]
  res <- cmb[, {
    .c <- enh_sequence
    for(i in seq(.N))
    {
      motWidth <- nchar(mot[i])
      is <- pos[i]-round(motWidth/2)
      ie <- is+motWidth-1
      substr(.c, is, ie) <- mot[i]
    }
    .c
  }, idx]$V1
  .(mutDref= unique(res))
}, .(BAID, group)]

#----------------------------#
# SAVE for Bernies' prediction
#----------------------------#
final <- rbindlist(list(WT= WT[, .(BAID, group, enh_sequence)],
                        add2Trl= addSynMot[, .(BAID, group, enh_sequence= addTrl)],
                        add3Trl= addSynMotInact[, .(BAID, group, enh_sequence= addTrl)],
                        add2Twist= addSynMot[, .(BAID, group, enh_sequence= addTwist)],
                        add3Twist= addSynMotInact[, .(BAID, group, enh_sequence= addTwist)],
                        add2Dref= addDrefMot[, .(BAID, group, enh_sequence= addDref)],
                        mutTrl= mutSynMot[group=="Trl", .(BAID, group, enh_sequence= mutSyn)],
                        mutTwist= mutSynMot[group=="twist", .(BAID, group, enh_sequence= mutSyn)],
                        mutDref= mutDrefMot[, .(BAID, group, enh_sequence= mutDref)]),
                   idcol= "mut")
final[mut=="WT", c("seqnames", "start", "end", "strand"):= .SD[WT, .(seqnames, start, end, strand), on="BAID"]]
final[, ID:= paste0(group, "__", mut, "__", rowid(group, mut, BAID), "__", BAID)]
seqinr::write.fasta(sequences = as.list(final$enh_sequence),
                    names = final$ID,
                    file.out = paste0("db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_with_WT.fasta"),
                    as.string = T)
fwrite(final, 
       "db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_with_WT.txt.gz",
       compress = "gzip")
