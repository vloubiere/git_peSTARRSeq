setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

dat <- fread("Rdata/metadata_processed.txt")
dat <- dat[DESeq2 & file.exists(pairs_counts)]
dat[as.data.table(read_excel("../../exp_data/vl_libraries.xlsx")), complexity:= i.Complexity, on= "vllib==lib_ID"]
dat <- dat[, fread(pairs_counts), .(vllib, complexity, pairs_counts, cdition)]
dat <- dat[, .(umi_counts= sum(umi_counts)), .(vllib, complexity, cdition, L, R)]
dat <- dat[, .(counts= sum(umi_counts), threshold= length(which(umi_counts>=5))), .(vllib, complexity, cdition)]

pdf("pdf/alignment/barplot_read_per_theoretical_pair.pdf")
par(mar= c(25,5,1,1))
barplot(dat$counts/dat$complexity, 
        col = sapply(dat$cdition, function(x) switch(x, 
                                                     "input"= "cornflowerblue",
                                                     "screen"= "tomato")),
        ylab= "reads/pair", 
        names.arg = paste0(dat$vllib,"_", dat$cdition), 
        las= 2,
        ylim= c(1, 2000),
        log= "y")
legend("topright",
       fill= c("cornflowerblue", "tomato"),
       legend= c("input", "screen_merged"),
       bty= "n")
barplot(dat$threshold/dat$complexity*100, 
        col = sapply(dat$cdition, function(x) switch(x, 
                                                     "input"= "cornflowerblue",
                                                     "screen"= "tomato")),
        ylab= "% pairs >= 5 reads", 
        names.arg = paste0(dat$vllib,"_", dat$cdition), 
        las= 2,
        ylim= c(0, 100))
legend("topright",
       fill= c("cornflowerblue", "tomato"),
       legend= c("input", "screen_merged"),
       bty= "n")
abline(h= 10)
dev.off()

file.show("pdf/alignment/barplot_read_per_theoretical_pair.pdf")