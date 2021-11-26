setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

if(!exists("dat"))
{
        dat <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
        dat <- as.data.table(dat)
        cols <- colnames(dat)
        dat[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]
        dat <- dat[, .(file= list.files("db/merged_counts/", 
                                        paste0(DESeq2_group, "_", cdition, "_rep", DESeq2_pseudo_rep, "_merged.txt"), full.names = T)),
                   .(DESeq2_group, cdition, DESeq2_pseudo_rep, Spike_in)]
        dat <- dat[, fread(file), (dat)]
        dat <- dat[type!="switched", .(umi_counts= sum(umi_counts)), .(cdition, DESeq2_group, L, R, Spike_in)]
}
res <- dat[, .(counts= sum(umi_counts), N_10_reads= .SD[umi_counts>=5, .N]), keyby= .(cdition, DESeq2_group, Spike_in)]
res[, lib:= tstrsplit(DESeq2_group, "_", keep= 1)]
res[as.data.table(read_excel("../../exp_data/vl_libraries.xlsx")), complexity:= i.Complexity, on= "lib==lib_ID"]
res[as.data.table(read_excel("../../exp_data/vl_libraries.xlsx")), complexity:= complexity+ifelse(is.na(i.Complexity), 0, i.Complexity), on= "Spike_in==lib_ID"]

pdf("pdf/alignment/barplot_read_per_theoretical_pair.pdf")
par(mar= c(25,5,1,1))
barplot(res$counts/res$complexity, 
        col = sapply(res$cdition, function(x) switch(x, 
                                                     "input"= "cornflowerblue",
                                                     "screen"= "tomato")),
        ylab= "reads/pair", 
        names.arg = res$DESeq2_group, 
        las= 2,
        ylim= c(1, 2000),
        log= "y")
legend("topright",
       fill= c("cornflowerblue", "tomato"),
       legend= c("input", "screen_merged"),
       bty= "n")
barplot(res$N_10_reads/res$complexity*100, 
        col = sapply(res$cdition, function(x) switch(x, 
                                                     "input"= "cornflowerblue",
                                                     "screen"= "tomato")),
        ylab= "% pairs >= 5 reads", 
        names.arg = res$DESeq2_group, 
        las= 2,
        ylim= c(0, 100))
legend("topright",
       fill= c("cornflowerblue", "tomato"),
       legend= c("input", "screen_merged"),
       bty= "n")
abline(h= 10)
dev.off()