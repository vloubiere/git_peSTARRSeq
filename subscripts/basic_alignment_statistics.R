setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)

meta <- read_excel("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
meta[, output_prefix:= paste0("/", my_ID, "__", gsub(".bam$", "", basename(BAM_path)))]

stats <- data.table(file= list.files("db/umi_counts/", "summary.txt$", full.names = T))
stats <- stats[, fread(file), file]
stats[, name:= gsub("_summary.txt", "", basename(file))]
stats[, Cc:= ifelse(meta[grep(name, output_prefix), DESeq2_group]=="NA", "red", "green"), name]
stats <- stats[order(grepl("screen", name), Cc, umi_collapsed_reads/total_reads)]

pdf("pdf/alignment/alignment_statistics.pdf", 
    height = nrow(stats)/5, 
    width = 10)
par(mar= c(7,30,2,2))
barplot(t(stats[, .(umi_collapsed_reads, total_reads)]), 
        beside = T,
        col= unlist(lapply(stats$Cc, function(x) c(x, "white"))), 
        names.arg = basename(stats$name),
        cex.names= 0.5, 
        # col.names= stats$Cc,
        horiz = T, 
        las= 1)
abline(v= 1e6, 
       lty= 2)
abline(v= 5e6, 
       lty= 2)
mtext(text = "N reads", 
      side = 1, 
      line = 5)
legend("bottomright", 
       bty= "n",
       fill= c("black", "white"), 
       legend = c("UMI collapsed", "all"))
dev.off()
