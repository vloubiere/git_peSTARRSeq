source("/groups/stark/vloubiere/scripts/R_functions/my_pheatmap_2.0.R")
require(data.table)

########### EXAMPLES ##############
HAM <- "dev_medium_C_00065"
SUP <- "dev_medium_C_00571"
SGL <- "dev_medium_B_00587"
ex <- data.table(enh_L= c("dev_strong_B_00023", "dev_strong_B_00289", "dev_medium_B_00522", "dev_medium_B_00591"), 
                 enh_R= c("dev_strong_C_00300", "dev_strong_C_00308", "dev_medium_C_00509", "dev_weak_C_00404"), 
                 ex= c("add", "add", "sup", "sup"))
motifs <- data.table(motif= c("flyfactorsurvey__srp_SANGER_5_FBgn0003507", "flyfactorsurvey__twi_da_SANGER_5_FBgn0000413", 
                              "jaspar__MA0491.1", "flyfactorsurvey__Trl_FlyReg_FBgn0013263", "homer__AVYTATCGATAD_DREF", 
                              "homer__MYGGTCACACTG_Unknown1"),
                     name= c("srp", "twist", "AP-1", "Trl", "DREF", "Ohler1"))

########### LOAD DATA ##############
if(!exists("ts") | !is.data.table(ts))
{
  ts <- fread("/groups/stark/vloubiere/projects/0002_library_design/Rdata/TWIST008_library_design_112019/template_switching_templates.txt")
}
if(!exists("lib"))
{
  lib <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/vl_library_112019.rds"))
}
if(!exists("feat"))
{
  feat <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/C_features_final_table.rds")
}
if(!exists("mot"))
{
  mot <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/D_motifs_final_table.rds")
}
if(!exists("counts_raw"))
{
  counts_raw <- readRDS("/groups/stark/vloubiere/data/pe_STARRSeq/Rdata/libvl002_SCR1_raw_counts.rds")
}
if(!exists("counts_clean"))
{
  counts_clean <- readRDS("/groups/stark/vloubiere/data/pe_STARRSeq/Rdata/libvl002_SCR1_clean_counts_table_final.rds")
}
if(!exists("dat"))
{
  dat <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/B_SCR1_peSTARRSeq_final_table.rds")
  dat <- dat[enh_L != enh_R]
  dat[median_L>2 & median_R>2 & grepl("^dev", enh_L) & grepl("^dev", enh_R), act_group:= "active~active"]
  dat[median_L>2 & median_R<0.25 & grepl("^dev", enh_L) & grepl("^dev", enh_R), act_group:= "active~inactive"]
  dat[median_L<0.25 & median_R>2 & grepl("^dev", enh_L) & grepl("^dev", enh_R), act_group:= "inactive~active"]
  dat[median_L<0.25 & median_R<0.25 & grepl("^dev", enh_L) & grepl("^dev", enh_R), act_group:= "inactive~inactive"]
  left_cuts <- unique(dat[median_L>=1 & grepl("^dev", enh_L) & grepl("^dev", enh_R), .(enh_L, median_L)])
  left_cuts <- c(-Inf, 0, quantile(left_cuts$median_L, seq(0, 1, length.out = 40)))
  dat[, quant_L:= cut(median_L, left_cuts, include.lowest = T)]
  right_cuts <- unique(dat[median_R>=1 & grepl("^dev", enh_L) & grepl("^dev", enh_R), .(enh_R, median_R)])
  right_cuts <- c(-Inf, 0, quantile(right_cuts$median_R, seq(0, 1, length.out = 40)))
  dat[, quant_R:= cut(median_R, right_cuts, include.lowest = T)]
  dat[diff > 1 & grepl("^dev", enh_L) & grepl("^dev", enh_R), diff_group:= "super-additive"]
  dat[diff > -1 & diff < 1 & grepl("^dev", enh_L) & grepl("^dev", enh_R), diff_group:= "additive"]
  dat[diff < -1 & grepl("^dev", enh_L) & grepl("^dev", enh_R), diff_group:= "sub-additive"]
}
if(!exists("ind"))
{
  ind <- unique(dat[!is.na(median_L), .(enh= enh_L, median_L)])
  ind <- merge(ind, unique(dat[!is.na(median_R), .(enh= enh_R, median_R)]), all.x= T, all.y= T)
}



