setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

meta <- data.table(file= list.files("db/umi_counts/", "vllib002", full.names = T))
meta[, c("lib", "cdition", "rep"):= tstrsplit(basename(file), "_|[.]", keep = 1:3)]
dat <- meta[, fread(file), .(cdition, rep)]
counts <- dat[, lapply(.SD, function(x) log2(sum(x)+1)), .(L, R, cdition), .SDcols= c("umi_counts", "total_counts")]

pdf("pdf/draft/umi_saturation.pdf", width = 3.5)
par(mfrow=c(2,1),
    mar= c(4,4,2,1),
    las= 1,
    tcl= -0.2,
    bty= "n",
    mgp= c(1.5,0.5,0))
counts[, {
  smoothScatter(total_counts,
                umi_counts, 
                xlab= "Total umi counts",
                ylab= "Collapsed umi counts",
                main= cdition)
  abline(0,1,lty= 2)
}, cdition]

par(mar= c(4,8,2,3),
    mgp= c(3.5,0.5,0))
umi_max <- dat[, max(umi_counts), .(cdition, rep)]
bar <- barplot(umi_max$V1,
        ylab= "Maximum number of UMIs",
        yaxt= "n")
ticks <- axTicks(2)
axis(2, 
     at= ticks, 
     labels = formatC(ticks, format= "e", digits = 0))
vl_tilt_xaxis(bar, 
              labels = umi_max[, paste0(cdition, " rep", rep)])
dev.off()

file.show("pdf/draft/umi_saturation.pdf")

