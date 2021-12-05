setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

stats <- fread("Rdata/metadata_processed.txt")
stats <- stats[file.exists(summary_counts) & DESeq2]
stats <- stats[, fread(summary_counts), .(my_ID, summary_counts)]
stats[, my_ID:= paste0(tstrsplit(my_ID, "_", keep= c(1, 9, 13)), collapse= "_"), my_ID]

pdf("pdf/alignment/aggregate_alignment_statistics.pdf", 
    height = nrow(stats)/3, 
    width = 15)
par(mar= c(7,15,2,2),
    mfrow= c(1, 2))
barplot(t(stats[, .(collapsed, mapped)]), 
        beside = T,
        col= c("green", "white"), 
        names.arg = basename(stats$my_ID),
        # col.names= stats$Cc,
        horiz = T, 
        las= 1, 
        log= "x")
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
pl <- t(stats[, .(collapsed/mapped*100)])
colnames(pl) <- t(stats[, my_ID])[1,]
barplot(pl,
        cex.names= 0.5,
        col= ifelse(grepl("input", colnames(pl)), "grey", "tomato"),
        beside = T,
        horiz = T, 
        las= 1, 
        xlim= c(0, 100),
        xlab= "% reads left after collapsing")
abline(v= 50, 
       lty= 2)
dev.off()
