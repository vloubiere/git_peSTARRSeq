setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(digest)

# Import
counts <- readRDS("Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2.rds")
counts <- counts[!grepl("temp_switch", rn)]
dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table_with_control_pairs.rds")
feat <- readRDS("Rdata/library/lib_features.rds")

#--------------------------#
# Apply cutoffs -> eligible pairs
#--------------------------#
sel <- data.table()

# Input counts cutoff
counts_cutoff_pdf <- "pdf/luciferase_validations/SECOND_ROUND_validations_input_counts_cutoff.pdf"
if(!file.exists(counts_cutoff_pdf))
{
  pdf(counts_cutoff_pdf)
  plot(NA, xlim= c(0, 8), ylim= c(0, 1))
  counts[, lapply(.SD, function(x) lines(density(log2(x+1)))), .SDcols= patterns("input_rep")]
  abline(v= 3.5, lty= 3)
  dev.off()
}

counts[, cutoff:= apply(.SD, 1, function(x) all(x>10)), .SDcols= patterns("input")]
sel <- dat[enh_L %in% counts[(cutoff), enh_L] & enh_R %in% counts[(cutoff), enh_R]]

# Select only dev enhancer and control pairs
sel <- na.omit(sel[enh_L %in% feat[group %in% c("dev", "control"), ID] & enh_R %in% feat[group %in% c("dev", "control"), ID]])

# Activity cutoffs
sel_ID <- sel[median_L>1.5 & median_L<3.5 & median_R>1.5 & median_R<3.5 & grepl("^dev", enh_L) & grepl("^dev", enh_R)]
# sel_ID <- sel_ID[enh_L %in% sel_ID[, any(diff>2) & any(diff<0.5), enh_L][(V1), enh_L] 
#                  & enh_R %in% sel_ID[, any(diff>2) & any(diff<0.5), enh_R][(V1), enh_R]]
sel_ctls_ID <- sel[(enh_L %in% sel_ID$enh_L & grepl("^control", enh_R)) | (grepl("^control", enh_L) & enh_R %in% sel_ID$enh_R)]
sel <- sel[enh_L %in% c(sel_ID$enh_L, sel_ctls_ID$enh_L) & enh_R %in% c(sel_ID$enh_R, sel_ctls_ID$enh_R)]

# # Select candidates with many pairs
sel <- sel[enh_L %in% sel[, .N, enh_L][N>20, enh_L]]
sel <- sel[enh_R %in% sel[, .N, enh_R][N>20, enh_R]]

# Select candidates with available homotypic pairs
sel <- sel[enh_L %in% enh_R]
sel <- sel[enh_R %in% enh_L]

# Print global heatmap
cpdf <- "pdf/luciferase_validations/SECOND_ROUND_heatmap_diff_for_candidaites_selection.pdf"
if(!file.exists(cpdf))
{
  sub <- sel[enh_L %in% sel[, .N, enh_L][N>150, enh_L]]
  sub <- sub[enh_R %in% sel[, .N, enh_R][N>100, enh_R]]
  sub <- sub[enh_L %in% enh_R]
  sub <- sub[enh_R %in% enh_L]
  
  pdf(cpdf, 20, 20)
  par(mar= c(10,10,10,10))
  my_heatmap(sub, "enh_L", "enh_R", "diff", cluster_rows = F, cluster_cols = F)
  dev.off()
}

#--------------------------#
# Compute many possible scenarios and select interesting ones
#--------------------------#
output <- "Rdata/luciferase_validations/SECOND_ROUND_training.rds"
if(!file.exists(output))
{
  cmb <- as.data.table(combinations(all, 4, replace=F))
  cmb[, idx:= .I]
  cmb <- melt(cmb, id.vars = "idx")
  cmb <- cmb[, CJ(enh_L= value, enh_R= value), idx]
  cmb <- sel[cmb, , on= c("enh_L", "enh_R"), nomatch= NA]
  cmb[, check:= .(!anyNA(diff) & any(diff<0.5)  & any(diff>2)), idx]
  cmb <- cmb[(check)]
  cmb[, c("median", "sd", "min", "max"):= .(median(diff), sd(diff), min(diff), max(diff)), idx]
  cmb <- cmb[order(-sd)]
  saveRDS(cmb, output)
}
cmb <- readRDS("Rdata/luciferase_validations/SECOND_ROUND_training.rds")

cpdf <- "pdf/luciferase_validations/SECOND_ROUND_final_sel_training.pdf"
if(!file.exists(cpdf))
{
  pdf(cpdf)
  cmb[idx %in% unique(cmb$idx)[1:500], {
    current <- my_heatmap(.SD, "enh_L", "enh_R", value.var = "diff", main = idx,
                          breaks = c(-2, -0.25, 0, 0.25, 4), col = c("cornflowerblue", "white", "white", "white", "tomato"))
    text(current$xcoor, current$ycoor, round(current$diff, 1))  
  }, idx]
  dev.off()
}

#--------------------------#
# FINAL SELECTION
#--------------------------#
final_sel <- c("dev_strong_B_00277", "dev_strong_B_00299", "dev_strong_C_00291", "dev_weak_C_00382", "dev_medium_C_00454", 
               "control_exon_C_00972", "control_exon_B_00960")
final <- unique(sel[enh_L %in% final_sel & enh_R %in% final_sel, .(enh_L, enh_R, diff)]) 
mat <- as.matrix(dcast.data.table(final, enh_L~enh_R, value.var = "diff"), 1)

# Heatmap
pdf("pdf/luciferase_validations/SECOND_ROUND_heatmap_SELECTED.pdf", width = 9, height = 8.5)
layout(matrix(c(2,1,4,3), ncol= 2, nrow= 2), widths = c(1,0.3), heights = c(0.25,1))

# Residuals
par(mar= c(12,12,2,4))
res <- my_heatmap(mat, cluster_rows = F, cluster_cols = F, breaks = c(-4,-0.5, 0, 0.5, 4), col = c("cornflowerblue", "white","white","white", "tomato"))
text(res$xcoor, res$ycoor, round(res$value, 1), cex= 0.6)

# Motifs Left
par(mar= c(0.5,12,2,4))
mot <- feat[ID %in% final$enh_R]
mot <- melt.data.table(mot, id.vars = c("ID", "motif_cl"), measure.vars = patterns("^motif__"))
mot <- mot[grepl("M0802$|GATA1$|AP1$|JUND_f1$", variable)]
mot[, value:= log2(value+1)]
mot <- mot[order(match(ID, res$variable))]
my_heatmap(mot, row.BY= "variable", col.BY= c("ID", "motif_cl"), value.var= "value", 
           cluster_cols = F, breaks = c(0,1,2), col = c("lightgrey", "white", "red"), col_labels = F)
# Motifs Right
par(mar= c(12,0.5,2,4))
mot <- feat[ID %in% final$enh_L]
mot <- melt.data.table(mot, id.vars = c("ID", "motif_cl"), measure.vars = patterns("^motif__"))
mot <- mot[grepl("M0802$|GATA1$|AP1$|JUND_f1$", variable)]
mot[, value:= log2(value+1)]
mot <- mot[order(match(ID, res$rn))]
my_heatmap(mot, row.BY = c("ID", "motif_cl"), col.BY = "variable", value.var = "value", 
           cluster_rows = F, breaks = c(0,1,2), col = c("lightgrey", "white", "red"), row_labels = F)
dev.off()

#------------------------#
# PRIMERS
#------------------------#
lib <- as.data.table(readRDS("Rdata/library/vl_library_112019.rds"))
lib <- lib[ID_vl %in% c(final$enh_L, final$enh_R)]
lib[, F_primer:= substr(enh_sequence, 1 , 22)]
lib[, R_primer:= as.character(reverseComplement(DNAStringSet(substr(enh_sequence, nchar(enh_sequence)-21 , nchar(enh_sequence)))))]
saveRDS(lib[, .(ID_vl, F_primer, R_primer)], "Rdata/luciferase_validations/SECOND_ROUND_final_candidates.rds")






