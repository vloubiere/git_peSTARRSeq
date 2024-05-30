setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Sup table 1
# Summary content libraries

# Sup table 2
dat <- readRDS("Rdata/vl_library_twist008_112019.rds")
dat <- as.data.table(dat)
dat <- dat[, .(ID= ID_vl, sub_library= linker_ID, seqnames, start, end, strand, genome,
               group, detail, fw_primer= fw_linker, enhancer_sequence= enh_sequence, rv_primer= rev_linker, oligo_full_sequence)]
tmp <- tempfile("WT_oligo_pool_sequences_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 4
dat <- readRDS("db/FC_tables/DSCP_large_WT_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("WT_oligo_pool_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 5
dat <- readRDS("db/FC_tables/DSCP_long_spacer_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("WT_long_spacer_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 6
dat <- readRDS("db/FC_tables/DSCP_ECD_WT_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("WT_ECD_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 7
dat <- readRDS("db/FC_tables/DSCP_OSC_WT_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("WT_OSC_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 8
mot <- readRDS("db/motif_counts/twist008_motif_counts_selected.rds")
mot <- vl_Dmel_motifs_DB_full[motif_ID %in% names(mot[, -1])]
mot <- mot[, {
  .c <- paste0("ID:", motif_ID, "|Cluster: ", motif_cluster, ifelse(is.na(Dmel), "", paste0("|TF:", Dmel)), "\n",
               paste0(paste0(c("A", "C", "G", "T"), ":", apply(as.matrix(mot[1, pwms_log_odds][[1]]), 1, paste0, collapse= ",")), collapse = "\n"))
}, motif_ID]
writeLines(paste0(mot$V1, collapse = "\n"),
           "Rdata/selected_motifs_sup_8.txt")
vl_dropbox_upload("Rdata/selected_motifs_sup_8.txt", ".")

# Sup table 9
dat <- readRDS("Rdata/vl_library_twist015_112022.rds")
dat <- as.data.table(dat)[sublib=="A"]
dat[, fw_primer:= "TTGACAGTGAGCGCGTCTCTCACCG"]
dat[, rv_primer:= "CCTAGGATCGACGCGGACAA"]
dat[, uniq_start:= seq(.N), gsub("^(.{110}).*", "\\1", enh_sequence)]
dat[, uniq_end:= seq(.N), gsub("^.*(.{110})$", "\\1", enh_sequence)]
dat[uniq_start==2, fw_primer:= paste0(fw_primer, "AT")]
dat[uniq_end==2, rv_primer:= paste0("AT", rv_primer)]
dat[, test:= paste0(fw_primer, enh_sequence, rv_primer)]
dat[nchar(test)==294, fw_primer:= paste0("GTA", fw_primer)]
dat[nchar(test)==294, rv_primer:= paste0(rv_primer, "ATG")]
dat[nchar(test)==296, fw_primer:= paste0("GT", fw_primer)]
dat[nchar(test)==296, rv_primer:= paste0(rv_primer, "TG")]
dat[nchar(test)==298, fw_primer:= paste0("G", fw_primer)]
dat[nchar(test)==298, rv_primer:= paste0(rv_primer, "G")]
dat[, test:= paste0(fw_primer, enh_sequence, rv_primer)]
dat[, identical(oligo_full_sequence, test)]
dat[, c("group", "detail", "seq_ID"):= tstrsplit(ID, "_", keep = c(1,2,4))]
dat <- dat[, .(ID, seqnames, start, end, strand, 
               group, detail, seq_ID, 
               fw_primer,
               enhancer_sequence= enh_sequence,
               rv_primer,
               oligo_full_sequence)]
tmp <- tempfile("Mutated_oligo_pool_sequences_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 10
dat <- readRDS("db/FC_tables/DSCP_mutant_library_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("Mutated_oligo_pool_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 11
dat <- readRDS("Rdata/vl_library_twist12_210610.rds")
dat <- as.data.table(dat)[linker_ID=="A"]
dat[, genome:= ifelse(detail=="ecoli", "Ecoli_20080805", "dm3")]
dat <- dat[, .(ID, seqnames, start, end, strand, genome,
               group, detail, fw_primer= fw_linker, enhancer_sequence= enh_seq, rv_primer= rev_linker, oligo_full_sequence)]
tmp <- tempfile("Focused_oligo_pool_sequences_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 12
dat <- readRDS("db/FC_tables/RpS12_focused_WT_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("Focused_oligo_pool_hkCP_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Sup table 13
dat <- readRDS("db/FC_tables/DSCP_focused_WT_FC_DESeq2.rds")
dat <- dat[, .(L, R, ctlL, ctlR, indL, padjL, indR, padjR, log2FoldChange, padj)]
tmp <- tempfile("Focused_oligo_pool_DSCP_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")
