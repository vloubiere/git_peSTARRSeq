# Alignment statistics ####
stats <- data.table(file= list.files("db/bam", "stats.txt", full.names = T))
stats <- stats[, .(Mapped_reads= readLines(file)), file]
stats <- stats[grepl("Uniquely_mapped_reads|Unmapped_reads", Mapped_reads)]
stats[, counts:= {current <- tstrsplit(Mapped_reads, " "); current[length(current)]}, Mapped_reads]
stats[, variable:= ifelse(grepl("Uniquely_mapped_reads", Mapped_reads), "Uniquely_mapped_reads", "Unmapped_reads"), Mapped_reads]
stats[, counts:= round(as.numeric(counts)/1e6, 2)]
stats <- dcast(stats, file~variable, value.var = "counts")
stats[, perc:= round(Uniquely_mapped_reads/(Uniquely_mapped_reads+Unmapped_reads)*100, 1)] 

pdf("pdf/alignment_statistics.pdf", width = 15, height = 30)
grid.table(stats)
dev.off()

