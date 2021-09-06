setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

stats <- data.table(file= list.files("db/merged_counts/", "summary.txt$", full.names = T))
stats <- stats[, fread(file), file]
stats[, name:= gsub("_merged.summary.txt", "", basename(file))]
stats <- stats[order(name, decreasing = T)]

pdf("pdf/alignment/aggregate_alignment_statistics.pdf", 
    height = nrow(stats)/3, 
    width = 15)
par(mar= c(7,15,2,2),
    mfrow= c(1, 2))
barplot(t(stats[, .(umi_collapsed, mapped)]), 
        beside = T,
        col= c("green", "white"), 
        names.arg = basename(stats$name),
        cex.names= 0.5, 
        # col.names= stats$Cc,
        horiz = T, 
        las= 1)
abline(v= 1e6, 
       lty= 2)
abline(v= 2.5e6, 
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
barplot(t(stats[, .(umi_collapsed/mapped*100)]),
        names.arg = basename(stats$name),
        cex.names= 0.5,
        col= stats[, ifelse(grepl("input", name), "grey", "tomato"), name]$V1,
        beside = T,
        horiz = T, 
        las= 1, 
        xlim= c(0, 100))
abline(v= 50, 
       lty= 2)
dev.off()
