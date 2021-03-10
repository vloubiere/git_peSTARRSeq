# Compute expected score ####
FC <- data.table(file= list.files("db/DE_analysis/", "DE.txt$", full.names = T))
FC[,  lib:= gsub("_DE.txt$", "", basename(file))]
FC[,  FC_file:= paste0("db/DE_analysis/", lib, "_add_scores.txt")]
FC[, {
  .c <- fread(file)
  .c[, spike_in:= ifelse(grepl("^ts", L) | grepl("^ts", R), T, F)]
  .c[, median_L:= .SD[grep("control", R), ifelse(.N<5, as.numeric(NA), median(log2FoldChange, na.rm = T))], L]
  .c[grepl("_C_", L) & R=="ts_SCR2_01002", median_L:= log2FoldChange, L]
  .c[, median_R:= .SD[grep("control", L), ifelse(.N<5, as.numeric(NA), median(log2FoldChange, na.rm = T))], R]
  .c[L=="ts_SCR2_01002" & grepl("_C_", R), median_R:= log2FoldChange, R]
  .c <- na.omit(.c[!grepl("control", L) & !grepl("control", R)])
  .c[, add:= log2(sum(2^median_L+2^median_R)), .(L, R)]
  fwrite(.c, FC_file, col.names = T, row.names = F, sep= "\t", quote= F)
}, .(file, FC_file)]

final <- FC[, fread(FC_file), FC]
final$file <- NULL
final$FC_file <- NULL
saveRDS(final, "Rdata/master_results_peSTARRSeq.rds")
####

