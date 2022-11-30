setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

if(!file.exists("db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_closest_activity_from_WT.rds"))
{
  # Import data
  pred <- fread("/groups/stark/almeida/Projects/Vincent_pairedSTARRseq/results/20221117_sequence_mutations_DeepSTARR_predictions/twist_0015_mutated_sequences_with_WT.fasta_DeepSTARR_predictions.txt")
  act <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
  dat <- fread("db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_with_WT.txt.gz")
  # Object
  dat[act, dev_log2FC:= dev_log2FoldChange, on= "BAID==ID"]
  dat[pred, dev_pred:= Predictions_dev, on= "ID==location"]
  dat[dat[mut=="WT"], WT_pred:= i.dev_pred, on= "BAID"]
  dat[, diff:= abs(dev_pred-WT_pred)]
  
  # Check bulk predicted vs observed activity
  check <- melt(dat[mut=="WT"], 
                id.vars = c("BAID", "group"),
                measure.vars = c("dev_log2FC", "dev_pred"))
  par(mar= c(12,4,2,2))
  boxplot(value~variable+group, check, las= 2, xlab= NA)
  
  # Select variants that are likely to have similar act
  dat <- dat[, .SD[which.min(diff)], .(mut, BAID)]
  dat$dev_log2FC <- dat$dev_pred <- dat$WT_pred <- NULL
  saveRDS(dat, 
          "db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_closest_activity_from_WT.rds")
}else
  dat <- readRDS("db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_closest_activity_from_WT.rds")

setorderv(dat, "diff")
dat[group=="control" & mut!="WT", .SD[1:50], mut]

# Select 40 ctls / motif pasting
ctlsAdd2Trl <- dat[group=="control" & mut=="add2Trl"][order(diff)][1:30]
ctlsAdd3Trl <- dat[group=="control" & mut=="add3Trl"][order(diff)][1:30]
ctlsAdd2Twist <- dat[group=="control" & mut=="add2Twist"][order(diff)][1:30]
ctlsAdd3Twist <- dat[group=="control" & mut=="add3Twist"][order(diff)][1:30]

# Select 40 DHS / motif pasting
DHSAdd2Trl <- dat[group=="DHS" & mut=="add2Trl"][order(diff)][1:30]
DHSAdd3Trl <- dat[group=="DHS" & mut=="add3Trl"][order(diff)][1:30]
DHSAdd2Twist <- dat[group=="DHS" & mut=="add2Twist"][order(diff)][1:30]
DHSAdd3Twist <- dat[group=="DHS" & mut=="add3Twist"][order(diff)][1:30]

# Select 50 noMotifAct noMotifAct / pasting motifs
noMotActAdd2Trl <- dat[group=="noMotifAct" & mut=="add2Trl"][order(diff)][1:50]
noMotActAdd2Twist <- dat[group=="noMotifAct" & mut=="add2Twist"][order(diff)][1:50]
noMotActAdd2Dref <- dat[group=="noMotifAct" & mut=="add2Dref"][order(diff)][1:35]

# Select 75 Trl enhancers that show low differences when mutating Trl
Trl <- dat[group=="Trl" & mut=="mutTrl"][order(diff)][1:50]
# Select 16 Trl enhancers that show low differences when pasting Dref
TrlAdd2Dref <- dat[BAID %in% Trl$BAID & mut=="add2Dref"]

# Select 75 Twist enhancers that show low differences when mutating Trl
Twist <- dat[group=="twist" & mut=="mutTwist"][1:50]
# Select 22 Twist enhancers that show low differences when pasting Dref
TwistAdd2Dref <- dat[BAID %in% Twist$BAID & mut=="add2Dref"]

# Select 50 shared enhancers that show low differences when mutating Dref
shared <- dat[group=="shared" & mut=="mutDref"][order(diff)][1:35]

#----------------------#
# Final object
#----------------------#
sel <- rbind(ctlsAdd2Trl,
             ctlsAdd3Trl,
             ctlsAdd2Twist,
             ctlsAdd3Twist,
             DHSAdd2Trl,
             DHSAdd3Trl,
             DHSAdd2Twist,
             DHSAdd3Twist,
             noMotActAdd2Trl,
             noMotActAdd2Twist,
             noMotActAdd2Dref,
             Trl,
             TrlAdd2Dref,
             Twist,
             TwistAdd2Dref,
             shared)
sel <- rbind(dat[mut=="WT" & BAID %in% sel$BAID],
             sel)
sel[, data.table(.N, .SD[, .N, keyby= mut]), keyby= group]

#-------------------------------------#
# Final objects and checks
#-------------------------------------#
lib <- sel[, .(ID)]
lib[, c("group", "mut", "BAID"):= tstrsplit(ID, "__", keep= c(1,2,4))]
lib[fread("db/STARR_predictions/twist_lib_015/twist_0015_mutated_sequences_with_WT.txt.gz"), 
    c("enh_sequence", "seqnames", "start", "end", "strand"):= .(i.enh_sequence, i.seqnames, i.start, i.end, i.strand), 
    on= "ID"]
# counts motifs
mot <- vl_motif_counts(lib$enh_sequence, 
                       sel= c("flyfactorsurvey__CG16778_SANGER_5_FBgn0003715", 
                              "cisbp__M2014",
                              "homer__AVYTATCGATAD_DREF"),
                       genome= "dm3", 
                       collapse_overlapping = T,
                       p.cutoff = 1e-4)
mot <- as.data.table(mot)
setnames(mot, c("motif_twist", "motif_Trl", "motif_Dref"))
lib <- cbind(lib, mot)
# Twist
lib[, .N, keyby= .(group, mut, motif_twist)]
vl_boxplot(list(controls= lib[group %in% c("control", "DHS", "noMotifAct") & mut=="WT", motif_twist],
                add2= lib[mut=="add2Twist", motif_twist], 
                add3= lib[mut=="add3Twist", motif_twist], 
                enr_cand= lib[group=="twist" & mut=="WT", motif_twist],
                mut= lib[mut=="mutTwist", motif_twist]), 
           outline= T)
# Trl
lib[, .N, keyby= .(group, mut, motif_Trl)]
vl_boxplot(list(controls= lib[group %in% c("control", "DHS", "noMotifAct") & mut=="WT", motif_Trl],
                add2= lib[mut=="add2Trl", motif_Trl], 
                add3= lib[mut=="add3Trl", motif_Trl], 
                enr_cand= lib[group=="Trl" & mut=="WT", motif_Trl],
                mut= lib[mut=="mutTrl", motif_Trl]), 
           outline= T)
lib[, motif_Trl:= sapply(stringr::str_locate_all(enh_sequence, "GAGAGAGA|CTCTCTCT"), nrow)]
vl_boxplot(list(controls= lib[group %in% c("control", "DHS", "noMotifAct") & mut=="WT", motif_Trl],
                add2= lib[mut=="add2Trl", motif_Trl], 
                add3= lib[mut=="add3Trl", motif_Trl], 
                enr_cand= lib[group=="Trl" & mut=="WT", motif_Trl],
                mut= lib[mut=="mutTrl", motif_Trl]), 
           outline= T)
# Dref
lib[, .N, keyby= .(group, mut, motif_Dref)]
vl_boxplot(list(controls= lib[group %in% c("noMotifAct", "Trl", "twist") & mut=="WT", motif_Dref],
                add2= lib[mut=="add2Dref", motif_Dref],
                enr_cand= lib[group=="shared" & mut=="WT", motif_Dref],
                mut= lib[mut=="mutDref", motif_Dref]), 
           outline= T)

saveRDS(lib, 
        "db/library_design/twist015/mutant_sublib_twist015.rds")
